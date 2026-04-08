import streamlit as st
import pandas as pd
import RNA
import peglit_min
import json
import base64
import asyncio
import signal
import time
from io import BytesIO
from contextlib import contextmanager

# ====================== 基础配置（适配5分钟+长序列） ======================
st.set_page_config(
    page_title="pegLIT - Long Sequence (5min+ Support)",
    layout="wide",
    initial_sidebar_state="expanded"
)

# 强制清除缓存（避免旧状态干扰）
st.cache_data.clear()
st.cache_resource.clear()

# ====================== 全局常量（适配长序列，支持5分钟+计算） ======================
DEFAULT_SEQ = {
    "spacer": "",
    "scaffold": "GTTTCAGAGCTATGCTGGAAACAGCATAGCAAGTTGAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGC",
    "template": "",
    "pbs": "",
    "linker": "NNNNNNNN",
    "motif": ""
}

# 长序列适配的长度限制（进一步放宽，满足更长序列需求）
SEQ_LENGTH_LIMITS = {
    "spacer": 100,    # 进一步放宽到100bp
    "template": 200,  # 进一步放宽到200bp
    "pbs": 80,        # 进一步放宽到80bp
    "motif": 500,     # 进一步放宽到500bp（支持更长motif）
    "linker": 100     # 进一步放宽到100bp
}

DEFAULT_PARAMS = {
    "ac_thresh": 0.5,
    "u_thresh": 3,
    "n_thresh": 3,
    "topn": 100,
    "epsilon": 1e-2,
    "num_repeats": 10,
    "num_steps": 250,
    "temp_init": 0.15,
    "temp_decay": 0.95,
    "bottleneck": 1,
    "seed": 2020
}

# 关键配置：长序列计算超时时间（适配5分钟以上，单位：秒）
BASE_TIMEOUT = 600  # 基础10分钟
MAX_TIMEOUT = 1800  # 最大30分钟
TIME_PER_100BP = 60  # 每100bp增加1分钟超时

# ====================== 会话状态初始化（新增计算耗时记录） ======================
if "rows" not in st.session_state:
    st.session_state.rows = [DEFAULT_SEQ.copy()]
if "params" not in st.session_state:
    st.session_state.params = DEFAULT_PARAMS.copy()
if "results" not in st.session_state:
    st.session_state.results = {}
if "calculated" not in st.session_state:
    st.session_state.calculated = False
if "calc_status" not in st.session_state:
    st.session_state.calc_status = "idle"  # idle/running/complete/error/timeout
if "calc_time" not in st.session_state:
    st.session_state.calc_time = {}  # 记录每行计算耗时

# ====================== 核心工具函数（重点优化超时和长序列计算） ======================
def update_seq_value(row_idx, seq_type):
    """实时更新序列值，添加长度过滤，避免超长输入"""
    if row_idx < len(st.session_state.rows):
        input_value = st.session_state.get(f"input_{seq_type}_{row_idx}", "").strip().upper()
        # 截断超长输入（避免前端限制失效，适配更长序列）
        max_len = SEQ_LENGTH_LIMITS.get(seq_type, 1000)
        st.session_state.rows[row_idx][seq_type] = input_value[:max_len]

@contextmanager
def timeout(seconds):
    """超时控制装饰器（适配5分钟+计算，兼容Windows/Linux）"""
    # Windows系统兼容方案（signal.SIGALRM不支持）
    if hasattr(signal, 'SIGALRM'):
        def timeout_handler(signum, frame):
            raise TimeoutError(f"Calculation exceeded {seconds//60}m {seconds%60}s")
        signal.signal(signal.SIGALRM, timeout_handler)
        signal.alarm(seconds)
        try:
            yield
        finally:
            signal.alarm(0)
    else:
        # Windows系统：用时间戳判断超时
        start_time = time.time()
        yield
        elapsed = time.time() - start_time
        if elapsed > seconds:
            raise TimeoutError(f"Calculation exceeded {seconds//60}m {seconds%60}s (actual: {elapsed//60:.0f}m {elapsed%60:.0f}s)")

async def calculate_single_row(row_data, params, row_idx, status_callback):
    """异步计算单行序列（重点适配5分钟+长序列，记录耗时）"""
    start_time = time.time()  # 记录开始时间
    try:
        # 基础值处理（保留原有逻辑，适配更长序列）
        spacer = row_data["spacer"].upper().strip()
        scaffold = row_data["scaffold"].upper().strip() or DEFAULT_SEQ["scaffold"]
        template = row_data["template"].upper().strip()
        pbs = row_data["pbs"].upper().strip()
        motif = row_data["motif"].upper().strip()
        linker_pattern = row_data["linker"].upper().strip() or DEFAULT_SEQ["linker"]
        
        # 更新状态，提示当前计算进度和序列长度
        seq_total_len = len(spacer) + len(template) + len(pbs) + len(motif)
        estimated_min = seq_total_len//100 + 5
        estimated_max = seq_total_len//100 + 10
        status_callback(f"Row {row_idx+1}: Total length {seq_total_len}bp | Estimated time: {estimated_min}~{estimated_max} minutes")
        
        # 动态计算超时时间（适配5分钟+，根据序列总长度调整）
        timeout_seconds = BASE_TIMEOUT + (seq_total_len // 100) * TIME_PER_100BP
        timeout_seconds = min(timeout_seconds, MAX_TIMEOUT)  # 不超过最大30分钟
        
        # 带超时的核心计算（适配长序列，允许5分钟以上耗时）
        with timeout(timeout_seconds):
            calc_result = peglit_min.pegLIT(
                seq_spacer=spacer,
                seq_scaffold=scaffold,
                seq_template=template,
                seq_pbs=pbs,
                seq_motif=motif,
                linker_pattern=linker_pattern,** params
            )
        
        # 解析结果（兼容不同返回格式）
        row_result = {"status": "success", "data": {}, "length_info": {}}
        best_linker = DEFAULT_SEQ["linker"]
        
        if isinstance(calc_result, dict):
            best_linker = calc_result.get("best_linker", calc_result.get("linker", best_linker))
            row_result["data"]["linker"] = best_linker
            row_result["data"]["candidates"] = calc_result.get("candidates", [])
        elif isinstance(calc_result, list) and len(calc_result) > 0:
            best_linker = calc_result[0] if isinstance(calc_result[0], str) else calc_result[0].get("linker", best_linker)
            row_result["data"]["linker"] = best_linker
        elif isinstance(calc_result, str):
            best_linker = calc_result
            row_result["data"]["linker"] = best_linker
        
        # 长序列二级结构预测优化（分长度处理，避免耗时过长）
        full_seq = f"{scaffold}{best_linker}{motif}"
        if len(full_seq) <= 500:  # 500bp以内正常预测
            ss, mfe = RNA.fold(full_seq)
            row_result["data"]["secondary_structure"] = ss
            row_result["data"]["mfe"] = mfe
        elif len(full_seq) <= 1000:  # 500-1000bp，简化预测（提速）
            RNA.cvar.uniq_ML = 1  # 启用简化模式
            ss, mfe = RNA.fold(full_seq)
            row_result["data"]["secondary_structure"] = f"{ss}\n(Simplified mode for long sequence)"
            row_result["data"]["mfe"] = mfe
        else:  # 1000bp以上，跳过预测
            row_result["data"]["secondary_structure"] = "Skipped (sequence > 1000bp, too long for RNAfold)"
            row_result["data"]["mfe"] = "N/A"
        
        # 记录长度信息和计算耗时
        end_time = time.time()
        calc_duration = round(end_time - start_time, 2)
        st.session_state.calc_time[row_idx] = calc_duration
        row_result["length_info"] = {
            "spacer": len(spacer),
            "template": len(template),
            "pbs": len(pbs),
            "motif": len(motif),
            "full_seq": len(full_seq),
            "calc_duration": calc_duration  # 计算耗时（秒）
        }
        
        return row_idx, row_result
    
    except TimeoutError as e:
        end_time = time.time()
        calc_duration = round(end_time - start_time, 2)
        st.session_state.calc_time[row_idx] = calc_duration
        return row_idx, {
            "status": "timeout", 
            "message": str(e), 
            "calc_duration": calc_duration
        }
    except Exception as e:
        end_time = time.time()
        calc_duration = round(end_time - start_time, 2)
        st.session_state.calc_time[row_idx] = calc_duration
        return row_idx, {
            "status": "error", 
            "message": f"Calculation failed: {str(e)}", 
            "calc_duration": calc_duration
        }

async def batch_calculate(rows, params, progress_callback, status_callback):
    """批量异步计算（支持多行长序列，每行长序列可耗时5分钟以上）"""
    tasks = []
    total_rows = len(rows)
    
    # 创建异步任务，每行独立计算，互不影响
    for idx, row in enumerate(rows):
        task = calculate_single_row(row, params, idx, status_callback)
        tasks.append(task)
    
    # 执行任务并更新进度
    results = {}
    for i, task in enumerate(asyncio.as_completed(tasks)):
        row_idx, row_result = await task
        results[row_idx] = row_result
        progress_callback((i+1)/total_rows)  # 实时更新进度
    
    return results

# ====================== 全局样式（新增长序列耗时/超时提示） ======================
st.markdown("""
<style>
/* 全局重置 */
* {font-family: 'Segoe UI', Roboto, Arial, sans-serif; box-sizing: border-box;}
#MainMenu, footer, header {visibility: hidden;}

/* 主容器 */
.main-container {max-width: 1400px; margin: 0 auto; padding: 20px;}

/* 标题区域 */
.title-section {
    text-align: center; margin-bottom: 30px; padding-bottom: 20px;
    border-bottom: 1px solid #e5e7eb;
}
.title-section h1 {font-size: 3.5rem; font-weight: 800; color: #111827; margin-bottom: 10px;}
.title-section p {font-size: 1.2rem; color: #4b5563; max-width: 800px; margin: 0 auto;}

/* 输入卡片 */
.input-card {
    background: #fff; border: 1px solid #e5e7eb; border-radius: 12px;
    padding: 24px; margin-bottom: 24px; box-shadow: 0 2px 4px rgba(0,0,0,0.05);
}

/* 表格样式 */
.input-table {width: 100%; border-collapse: collapse; margin: 16px 0;}
.input-table th {
    padding: 12px 16px; text-align: left; font-size: 1.1rem;
    font-weight: 600; color: #1f2937; border-bottom: 2px solid #e5e7eb;
}
.input-table td {
    padding: 12px 16px; border-bottom: 1px solid #e5e7eb;
    vertical-align: top;
}
.input-table input {
    width: 100%; padding: 8px 12px; border: 1px solid #d1d5db;
    border-radius: 6px; font-size: 1rem; outline: none; margin-top: 6px;
}
.input-table input:focus {
    border-color: #2563eb; box-shadow: 0 0 0 2px rgba(37, 99, 235, 0.1);
}

/* 标签样式 */
.seq-label {font-size: 0.9rem; font-weight: 500; color: #4b5563; display: block;}
.required {color: #dc2626; margin-left: 4px;}
.length-hint {font-size: 0.8rem; color: #6b7280; margin-left: 4px;}

/* 按钮样式 */
.stButton > button {border-radius: 8px !important;}
.stButton > button[data-testid="baseButton-primary"] {
    background-color: #2563eb !important; color: white !important;
    padding: 12px 24px !important; font-size: 1.1rem !important;
    font-weight: 600 !important;
}
.stButton > button[data-testid="baseButton-secondary"] {
    background-color: #f3f4f6 !important; color: #1f2937 !important;
    border: 1px solid #d1d5db !important; padding: 8px 16px !important;
    font-size: 0.9rem !important;
}

/* 进度提示 */
.progress-container {
    margin: 20px 0; padding: 16px; border-radius: 8px;
    background-color: #f9fafb; border: 1px solid #e5e7eb;
}
.progress-container strong {color: #2563eb;}

/* 结果卡片 */
.result-card {
    background: #fff; border: 1px solid #e5e7eb; border-radius: 12px;
    padding: 24px; margin-top: 24px; box-shadow: 0 2px 4px rgba(0,0,0,0.05);
}
.result-card h3 {
    font-size: 1.5rem; font-weight: 600; color: #1f2937;
    margin-bottom: 16px; padding-bottom: 8px; border-bottom: 1px solid #e5e7eb;
}
.result-card .length-info {font-size: 0.9rem; color: #6b7280; margin-bottom: 8px;}
.result-card .calc-time {font-size: 0.9rem; color: #2563eb; margin-bottom: 12px;}

/* 二级结构展示 */
.rna-structure {
    font-family: monospace; font-size: 0.9rem; line-height: 1.6;
    padding: 16px; background-color: #f9fafb; border-radius: 8px;
    white-space: pre-wrap; margin: 16px 0; max-height: 400px;
    overflow-y: auto; border: 1px solid #e5e7eb;
}

/* 长序列提示框 */
.long-seq-hint {
    background-color: #fffbeb; border-left: 4px solid #f59e0b;
    padding: 12px 16px; margin: 16px 0; border-radius: 4px;
}

/* 超时/错误提示 */
.timeout-alert {background-color: #fff2f2; border-left: 4px solid #f87171; padding: 12px; border-radius: 4px; margin: 8px 0;}
.error-alert {background-color: #fef2f2; border-left: 4px solid #dc2626; padding: 12px; border-radius: 4px; margin: 8px 0;}
</style>
""", unsafe_allow_html=True)

# ====================== 侧边栏参数配置 ======================
with st.sidebar:
    st.markdown("<h2 style='font-size:1.8rem; font-weight:700; margin-bottom:20px;'>Parameters</h2>", unsafe_allow_html=True)
    
    # 能量阈值参数
    st.session_state.params["ac_thresh"] = st.number_input(
        "AC Threshold", min_value=0.0, max_value=1.0,
        value=DEFAULT_PARAMS["ac_thresh"], step=0.05,
        help="Acceptable threshold for AC content (lower = stricter)"
    )
    
    st.session_state.params["u_thresh"] = st.number_input(
        "U Threshold", min_value=0, max_value=10,
        value=DEFAULT_PARAMS["u_thresh"], step=1,
        help="Maximum allowed consecutive U's (higher = more permissive)"
    )
    
    st.session_state.params["n_thresh"] = st.number_input(
        "N Threshold", min_value=0, max_value=10,
        value=DEFAULT_PARAMS["n_thresh"], step=1,
        help="Maximum allowed consecutive N's (higher = more permissive)"
    )
    
    # 模拟退火参数（长序列可减少迭代次数提速）
    st.markdown("<h3 style='font-size:1.2rem; font-weight:600; margin:20px 0 10px;'>Simulated Annealing</h3>", unsafe_allow_html=True)
    
    st.session_state.params["topn"] = st.number_input(
        "Top N (candidates)", min_value=10, max_value=500,
        value=DEFAULT_PARAMS["topn"], step=10,
        help="Reduce for long sequences (faster calculation)"
    )
    
    st.session_state.params["num_steps"] = st.number_input(
        "Number of Steps", min_value=50, max_value=1000,
        value=DEFAULT_PARAMS["num_steps"], step=50,
        help="Reduce for long sequences (faster calculation)"
    )
    
    st.session_state.params["num_repeats"] = st.number_input(
        "Number of Repeats", min_value=1, max_value=50,
        value=DEFAULT_PARAMS["num_repeats"], step=1,
        help="Reduce for long sequences (faster calculation)"
    )
    
    # 随机种子
    st.session_state.params["seed"] = st.number_input(
        "Random Seed", min_value=1, max_value=9999,
        value=DEFAULT_PARAMS["seed"], step=1,
        help="Random seed for reproducibility"
    )

# ====================== 主界面 ======================
st.markdown("<div class='main-container'>", unsafe_allow_html=True)

# 标题区域
st.markdown("""
<div class="title-section">
    <h1>pegLIT (Long Sequence Support)</h1>
    <p>Automatically identify non-interfering nucleotide linkers (supports 5min+ long sequence calculation)</p>
</div>
""", unsafe_allow_html=True)

# 长序列提示
st.markdown("""
<div class="long-seq-hint">
    <strong>💡 Long Sequence Tips:</strong> 
    1. Calculation time for long sequences (>200bp) may exceed 5 minutes; 
    2. RNAfold is simplified/skipped for sequences >500bp; 
    3. Reduce 'Top N/Number of Steps' in sidebar for faster results;
    4. Maximum calculation time per row: 30 minutes.
</div>
""", unsafe_allow_html=True)

# 输入区域
st.markdown("<div class='input-card'>", unsafe_allow_html=True)
st.markdown("<h3 style='font-size:1.5rem; font-weight:600; margin-bottom:16px;'>Sequence Input</h3>", unsafe_allow_html=True)

# 渲染输入表格
st.markdown("<table class='input-table'>", unsafe_allow_html=True)
st.markdown("""
<thead>
    <tr>
        <th>Spacer <span class='required'>*</span> <span class='length-hint'>(≤100bp)</span></th>
        <th>Scaffold</th>
        <th>Template <span class='required'>*</span> <span class='length-hint'>(≤200bp)</span></th>
        <th>PBS <span class='required'>*</span> <span class='length-hint'>(≤80bp)</span></th>
        <th>Linker Pattern <span class='length-hint'>(≤100bp)</span></th>
        <th>Motif <span class='required'>*</span> <span class='length-hint'>(≤500bp)</span></th>
    </tr>
</thead>
<tbody>
""", unsafe_allow_html=True)

# 遍历渲染每一行输入框（带长度限制）
for row_idx, row_data in enumerate(st.session_state.rows):
    st.markdown("<tr>", unsafe_allow_html=True)
    
    # 1. Spacer列（必填，长序列放宽到100bp）
    st.markdown("<td>", unsafe_allow_html=True)
    st.markdown("<span class='seq-label'>Spacer <span class='required'>*</span></span>", unsafe_allow_html=True)
    st.text_input(
        label=f"spacer_{row_idx}",
        value=row_data["spacer"],
        label_visibility="collapsed",
        key=f"input_spacer_{row_idx}",
        max_chars=SEQ_LENGTH_LIMITS["spacer"],
        on_change=update_seq_value,
        args=(row_idx, "spacer")
    )
    st.markdown("</td>", unsafe_allow_html=True)
    
    # 2. Scaffold列（可选）
    st.markdown("<td>", unsafe_allow_html=True)
    st.markdown("<span class='seq-label'>Scaffold</span>", unsafe_allow_html=True)
    st.text_input(
        label=f"scaffold_{row_idx}",
        value=row_data["scaffold"],
        label_visibility="collapsed",
        key=f"input_scaffold_{row_idx}",
        on_change=update_seq_value,
        args=(row_idx, "scaffold")
    )
    st.markdown("</td>", unsafe_allow_html=True)
    
    # 3. Template列（必填，长序列放宽到200bp）
    st.markdown("<td>", unsafe_allow_html=True)
    st.markdown("<span class='seq-label'>Template <span class='required'>*</span></span>", unsafe_allow_html=True)
    st.text_input(
        label=f"template_{row_idx}",
        value=row_data["template"],
        label_visibility="collapsed",
        key=f"input_template_{row_idx}",
        max_chars=SEQ_LENGTH_LIMITS["template"],
        on_change=update_seq_value,
        args=(row_idx, "template")
    )
    st.markdown("</td>", unsafe_allow_html=True)
    
    # 4. PBS列（必填，长序列放宽到80bp）
    st.markdown("<td>", unsafe_allow_html=True)
    st.markdown("<span class='seq-label'>PBS <span class='required'>*</span></span>", unsafe_allow_html=True)
    st.text_input(
        label=f"pbs_{row_idx}",
        value=row_data["pbs"],
        label_visibility="collapsed",
        key=f"input_pbs_{row_idx}",
        max_chars=SEQ_LENGTH_LIMITS["pbs"],
        on_change=update_seq_value,
        args=(row_idx, "pbs")
    )
    st.markdown("</td>", unsafe_allow_html=True)
    
    # 5. Linker Pattern列（可选，长序列放宽到100bp）
    st.markdown("<td>", unsafe_allow_html=True)
    st.markdown("<span class='seq-label'>Linker Pattern</span>", unsafe_allow_html=True)
    st.text_input(
        label=f"linker_{row_idx}",
        value=row_data["linker"],
        label_visibility="collapsed",
        disabled=st.session_state.calculated,
        key=f"input_linker_{row_idx}",
        max_chars=SEQ_LENGTH_LIMITS["linker"],
        on_change=update_seq_value,
        args=(row_idx, "linker")
    )
    st.markdown("</td>", unsafe_allow_html=True)
    
    # 6. Motif列（必填，长序列放宽到500bp）
    st.markdown("<td>", unsafe_allow_html=True)
    st.markdown("<span class='seq-label'>Motif <span class='required'>*</span></span>", unsafe_allow_html=True)
    st.text_input(
        label=f"motif_{row_idx}",
        value=row_data["motif"],
        label_visibility="collapsed",
        key=f"input_motif_{row_idx}",
        max_chars=SEQ_LENGTH_LIMITS["motif"],
        on_change=update_seq_value,
        args=(row_idx, "motif")
    )
    st.markdown("</td>", unsafe_allow_html=True)
    
    st.markdown("</tr>", unsafe_allow_html=True)

st.markdown("</tbody></table>", unsafe_allow_html=True)

# 功能按钮行（Add/Import/Export）
st.markdown("<div class='func-btn-row'>", unsafe_allow_html=True)
col_add, col_import, col_export = st.columns([1, 1, 1])
with col_add:
    if st.button("➕ Add New Row", key="add_row", type="secondary", help="Add new sequence row"):
        st.session_state.rows.append(DEFAULT_SEQ.copy())
        st.rerun()

with col_import:
    uploaded_file = st.file_uploader(
        "Import CSV", type="csv", label_visibility="collapsed", key="csv_upload"
    )
    if uploaded_file is not None:
        try:
            df = pd.read_csv(uploaded_file)
            df.columns = df.columns.str.lower()
            required_cols = ["spacer", "scaffold", "template", "pbs", "linker", "motif"]
            for col in required_cols:
                if col not in df.columns:
                    df[col] = DEFAULT_SEQ[col]
            df = df[required_cols]
            
            # 截断超长序列（兼容导入的长序列）
            for col in SEQ_LENGTH_LIMITS.keys():
                if col in df.columns:
                    df[col] = df[col].astype(str).str.strip().str.upper().str[:SEQ_LENGTH_LIMITS[col]]
            
            st.session_state.rows = []
            for _, row in df.iterrows():
                st.session_state.rows.append(row.to_dict())
            st.success("✅ CSV imported successfully (long sequences truncated to limits)")
            st.rerun()
        except Exception as e:
            st.error(f"❌ CSV import failed: {str(e)}")

with col_export:
    def export_to_csv():
        df = pd.DataFrame(st.session_state.rows)
        buffer = BytesIO()
        df.to_csv(buffer, index=False, encoding="utf-8")
        buffer.seek(0)
        return base64.b64encode(buffer.getvalue()).decode()
    
    if st.session_state.rows:
        csv_base64 = export_to_csv()
        st.download_button(
            label="📤 Export CSV",
            data=base64.b64decode(csv_base64),
            file_name="peglit_long_sequences.csv",
            mime="text/csv",
            type="secondary",
            key="export_csv"
        )
st.markdown("</div>", unsafe_allow_html=True)

# 计算按钮 + 异步计算逻辑（核心适配长序列）
st.markdown("<div style='text-align:center; margin:24px 0;'>", unsafe_allow_html=True)
if st.button("START CALCULATION (Long Seq Support)", key="start_calc", type="primary", 
            disabled=st.session_state.calc_status == "running"):
    
    # 前置校验（必填+长度）
    validation_error = False
    error_msg = "❌ Please fix the following issues:\n"
    
    for row_idx, row_data in enumerate(st.session_state.rows):
        missing_fields = []
        length_warnings = []
        
        # 检查必填字段
        if not row_data["spacer"].strip():
            missing_fields.append("Spacer")
        if not row_data["template"].strip():
            missing_fields.append("Template")
        if not row_data["pbs"].strip():
            missing_fields.append("PBS")
        if not row_data["motif"].strip():
            missing_fields.append("Motif")
        
        # 检查长度（仅警告，不阻止计算）
        for seq_type, max_len in SEQ_LENGTH_LIMITS.items():
            current_len = len(row_data[seq_type].strip())
            if current_len > max_len:
                length_warnings.append(f"{seq_type} (truncated to {max_len}bp)")
        
        # 收集错误
        if missing_fields:
            validation_error = True
            error_msg += f"- Row {row_idx+1}: Missing fields - {', '.join(missing_fields)}\n"
        if length_warnings:
            st.warning(f"⚠️ Row {row_idx+1}: {', '.join(length_warnings)}")
    
    if validation_error:
        st.error(error_msg)
    else:
        # 初始化计算状态
        st.session_state.calc_status = "running"
        st.session_state.results = {}
        st.session_state.calculated = False
        
        # 创建进度和状态展示
        progress_container = st.markdown("<div class='progress-container'><strong>Status:</strong> Initializing calculation...</div>", unsafe_allow_html=True)
        progress_bar = st.progress(0)
        status_text = st.empty()
        
        # 异步执行计算
        try:
            # 适配Windows系统的asyncio事件循环
            try:
                loop = asyncio.get_event_loop()
            except RuntimeError:
                loop = asyncio.new_event_loop()
                asyncio.set_event_loop(loop)
            
            # 定义进度/状态回调
            def update_progress(progress):
                progress_bar.progress(progress)
            
            def update_status(msg):
                status_text.text(msg)
                progress_container.markdown(f"<div class='progress-container'><strong>Status:</strong> {msg}</div>", unsafe_allow_html=True)
            
            # 执行批量计算
            results = loop.run_until_complete(
                batch_calculate(
                    st.session_state.rows,
                    st.session_state.params,
                    update_progress,
                    update_status
                )
            )
            
            # 保存结果
            st.session_state.results = results
            st.session_state.calculated = True
            st.session_state.calc_status = "complete"
            
            # 清理进度
            progress_bar.empty()
            status_text.empty()
            progress_container.empty()
            
            st.success("✅ Calculation completed (long sequence mode)!")
            st.rerun()
            
        except Exception as e:
            st.error(f"❌ Calculation failed: {str(e)}")
            st.session_state.calc_status = "error"
            progress_bar.empty()
            status_text.empty()
            progress_container.empty()
        finally:
            try:
                loop.close()
            except:
                pass

st.markdown("</div>", unsafe_allow_html=True)
st.markdown("</div>", unsafe_allow_html=True)  # 关闭input-card

# 结果展示区域（适配长序列）
if st.session_state.calculated and st.session_state.results:
    st.markdown("<div class='result-card'>", unsafe_allow_html=True)
    st.markdown("<h3>Calculation Results (Long Sequence)</h3>", unsafe_allow_html=True)
    
    total_calc_time = sum(st.session_state.calc_time.values())
    st.markdown(f"<div style='font-size:0.95rem; color:#4b5563; margin-bottom:20px;'>"
                f"<strong>Total calculation time:</strong> {total_calc_time//60:.0f}m {total_calc_time%60:.0f}s "
                f"({len(st.session_state.rows)} rows processed)</div>", unsafe_allow_html=True)
    
    for row_idx in sorted(st.session_state.results.keys()):
        result = st.session_state.results[row_idx]
        st.markdown(f"<h4 style='margin:20px 0 10px;'>Row {row_idx+1}</h4>", unsafe_allow_html=True)
        
        # 计算耗时展示
        calc_duration = result.get("calc_duration", 0)
        st.markdown(f"""
        <div class='calc-time'>
            <strong>Calculation Time:</strong> {calc_duration//60:.0f}m {calc_duration%60:.0f}s
        </div>
        """, unsafe_allow_html=True)
        
        # 超时结果
        if result["status"] == "timeout":
            st.markdown(f"""
            <div class='timeout-alert'>
                <strong>⚠️ Timeout Error:</strong> {result.get('message', 'Unknown timeout')}
                <br><small>Tip: Increase 'BASE_TIMEOUT' in code or reduce sequence length</small>
            </div>
            """, unsafe_allow_html=True)
            continue
        
        # 错误结果
        if result["status"] == "error":
            st.markdown(f"""
            <div class='error-alert'>
                <strong>❌ Calculation Error:</strong> {result.get('message', 'Unknown error')}
            </div>
            """, unsafe_allow_html=True)
            continue
        
        # 成功结果
        result_data = result["data"]
        length_info = result.get("length_info", {})
        
        # 长度信息展示
        st.markdown(f"""
        <div class='length-info'>
            <strong>Sequence Lengths:</strong> 
            Spacer={length_info.get('spacer', 0)}bp | 
            Template={length_info.get('template', 0)}bp | 
            PBS={length_info.get('pbs', 0)}bp | 
            Motif={length_info.get('motif', 0)}bp | 
            Full Seq={length_info.get('full_seq', 0)}bp
        </div>
        """, unsafe_allow_html=True)
        
        # 核心结果
        col1, col2 = st.columns([1, 2])
        with col1:
            st.markdown("<strong>Optimal Linker:</strong>", unsafe_allow_html=True)
            st.markdown(f"<div style='font-size:1.2rem; font-weight:600; color:#2563eb;'>{result_data['linker']}</div>", unsafe_allow_html=True)
            
            st.markdown("<strong>Minimum Free Energy (MFE):</strong>", unsafe_allow_html=True)
            st.markdown(f"<div style='font-size:1rem;'>{result_data['mfe']}</div>", unsafe_allow_html=True)
        
        with col2:
            st.markdown("<strong>RNA Secondary Structure:</strong>", unsafe_allow_html=True)
            st.markdown(f"<div class='rna-structure'>{result_data['secondary_structure']}</div>", unsafe_allow_html=True)
        
        # 候选Linker（长序列只展示前10个，避免渲染卡顿）
        if "candidates" in result_data and result_data["candidates"]:
            st.markdown("<strong>Top Candidate Linkers (Filtered for Long Seq):</strong>", unsafe_allow_html=True)
            candidates_df = pd.DataFrame(result_data["candidates"])
            st.dataframe(
                candidates_df.head(10),
                use_container_width=True,
                hide_index=True
            )
        
        st.markdown("<hr style='margin:20px 0; border:1px solid #e5e7eb;'>", unsafe_allow_html=True)
    
    st.markdown("</div>", unsafe_allow_html=True)

st.markdown("</div>", unsafe_allow_html=True)  # 关闭main-container
