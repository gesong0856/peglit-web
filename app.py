import streamlit as st
import pandas as pd
import RNA
import peglit_min
import json
import base64
from io import BytesIO

# ====================== 基础配置 ======================
st.set_page_config(
    page_title="pegLIT - Official",
    layout="wide",
    initial_sidebar_state="expanded"
)

# 强制清除缓存
st.cache_data.clear()
st.cache_resource.clear()

# ====================== 全局常量 ======================
DEFAULT_SEQ = {
    "spacer": "",
    "scaffold": "GTTTCAGAGCTATGCTGGAAACAGCATAGCAAGTTGAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGC",
    "template": "",
    "pbs": "",
    "linker": "NNNNNNNN",
    "motif": ""
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

# ====================== 会话状态初始化 ======================
if "rows" not in st.session_state:
    st.session_state.rows = [DEFAULT_SEQ.copy()]
if "params" not in st.session_state:
    st.session_state.params = DEFAULT_PARAMS.copy()
if "results" not in st.session_state:
    st.session_state.results = {}
if "calculated" not in st.session_state:
    st.session_state.calculated = False

# ====================== 全局样式（复刻官网） ======================
st.markdown("""
<style>
/* 全局重置 */
* {
    font-family: 'Segoe UI', Roboto, Arial, sans-serif;
    box-sizing: border-box;
}

/* 隐藏默认元素 */
#MainMenu, footer, header {visibility: hidden;}

/* 主容器 */
.main-container {
    max-width: 1400px;
    margin: 0 auto;
    padding: 20px;
}

/* 标题区域 */
.title-section {
    text-align: center;
    margin-bottom: 30px;
    padding-bottom: 20px;
    border-bottom: 1px solid #e5e7eb;
}
.title-section h1 {
    font-size: 3.5rem;
    font-weight: 800;
    color: #111827;
    margin-bottom: 10px;
}
.title-section p {
    font-size: 1.2rem;
    color: #4b5563;
    max-width: 800px;
    margin: 0 auto;
}

/* 输入卡片 */
.input-card {
    background: #ffffff;
    border: 1px solid #e5e7eb;
    border-radius: 12px;
    padding: 24px;
    margin-bottom: 24px;
    box-shadow: 0 2px 4px rgba(0,0,0,0.05);
}

/* 表格样式 */
.input-table {
    width: 100%;
    border-collapse: collapse;
    margin: 16px 0;
}
.input-table th {
    padding: 12px 16px;
    text-align: left;
    font-size: 1.1rem;
    font-weight: 600;
    color: #1f2937;
    border-bottom: 2px solid #e5e7eb;
}
.input-table td {
    padding: 12px 16px;
    border-bottom: 1px solid #e5e7eb;
}
.input-table input {
    width: 100%;
    padding: 8px 12px;
    border: 1px solid #d1d5db;
    border-radius: 6px;
    font-size: 1rem;
    outline: none;
}
.input-table input:focus {
    border-color: #2563eb;
    box-shadow: 0 0 0 2px rgba(37, 99, 235, 0.1);
}
.input-table input:disabled {
    background-color: #f9fafb;
    color: #4b5563;
    cursor: not-allowed;
}

/* 按钮样式 */
.btn-primary {
    background-color: #2563eb !important;
    color: white !important;
    border: none !important;
    border-radius: 8px !important;
    padding: 12px 24px !important;
    font-size: 1.1rem !important;
    font-weight: 600 !important;
    cursor: pointer !important;
    transition: background-color 0.2s !important;
}
.btn-primary:hover {
    background-color: #1d4ed8 !important;
}
.btn-secondary {
    background-color: #f3f4f6 !important;
    color: #1f2937 !important;
    border: 1px solid #d1d5db !important;
    border-radius: 8px !important;
    padding: 8px 16px !important;
    font-size: 0.9rem !important;
    cursor: pointer !important;
}
.btn-secondary:hover {
    background-color: #e5e7eb !important;
}

/* 结果卡片 */
.result-card {
    background: #ffffff;
    border: 1px solid #e5e7eb;
    border-radius: 12px;
    padding: 24px;
    margin-top: 24px;
    box-shadow: 0 2px 4px rgba(0,0,0,0.05);
}
.result-card h3 {
    font-size: 1.5rem;
    font-weight: 600;
    color: #1f2937;
    margin-bottom: 16px;
    padding-bottom: 8px;
    border-bottom: 1px solid #e5e7eb;
}

/* 二级结构展示 */
.rna-structure {
    font-family: monospace;
    font-size: 1.1rem;
    line-height: 1.6;
    padding: 16px;
    background-color: #f9fafb;
    border-radius: 8px;
    white-space: pre-wrap;
    margin: 16px 0;
}

/* 侧边栏样式 */
.sidebar-param {
    margin-bottom: 16px;
}
.sidebar-param label {
    font-size: 0.95rem;
    font-weight: 500;
    color: #1f2937;
    margin-bottom: 4px;
    display: block;
}
.sidebar-param input {
    width: 100%;
    padding: 8px;
    border-radius: 6px;
    border: 1px solid #d1d5db;
}
</style>
""", unsafe_allow_html=True)

# ====================== 侧边栏参数配置（复刻官网） ======================
with st.sidebar:
    st.markdown("<h2 style='font-size:1.8rem; font-weight:700; margin-bottom:20px;'>Parameters</h2>", unsafe_allow_html=True)
    
    # 能量阈值参数
    st.markdown("<div class='sidebar-param'>", unsafe_allow_html=True)
    st.session_state.params["ac_thresh"] = st.number_input(
        "AC Threshold",
        min_value=0.0,
        max_value=1.0,
        value=DEFAULT_PARAMS["ac_thresh"],
        step=0.05,
        help="Acceptable threshold for AC content"
    )
    st.markdown("</div>", unsafe_allow_html=True)
    
    st.markdown("<div class='sidebar-param'>", unsafe_allow_html=True)
    st.session_state.params["u_thresh"] = st.number_input(
        "U Threshold",
        min_value=0,
        max_value=10,
        value=DEFAULT_PARAMS["u_thresh"],
        step=1,
        help="Maximum allowed consecutive U's"
    )
    st.markdown("</div>", unsafe_allow_html=True)
    
    st.markdown("<div class='sidebar-param'>", unsafe_allow_html=True)
    st.session_state.params["n_thresh"] = st.number_input(
        "N Threshold",
        min_value=0,
        max_value=10,
        value=DEFAULT_PARAMS["n_thresh"],
        step=1,
        help="Maximum allowed consecutive N's"
    )
    st.markdown("</div>", unsafe_allow_html=True)
    
    # 模拟退火参数
    st.markdown("<h3 style='font-size:1.2rem; font-weight:600; margin:20px 0 10px;'>Simulated Annealing</h3>", unsafe_allow_html=True)
    
    st.markdown("<div class='sidebar-param'>", unsafe_allow_html=True)
    st.session_state.params["topn"] = st.number_input(
        "Top N",
        min_value=10,
        max_value=500,
        value=DEFAULT_PARAMS["topn"],
        step=10,
        help="Number of top candidates to keep"
    )
    st.markdown("</div>", unsafe_allow_html=True)
    
    st.markdown("<div class='sidebar-param'>", unsafe_allow_html=True)
    st.session_state.params["temp_init"] = st.number_input(
        "Initial Temperature",
        min_value=0.05,
        max_value=0.5,
        value=DEFAULT_PARAMS["temp_init"],
        step=0.01,
        help="Initial temperature for annealing"
    )
    st.markdown("</div>", unsafe_allow_html=True)
    
    st.markdown("<div class='sidebar-param'>", unsafe_allow_html=True)
    st.session_state.params["temp_decay"] = st.number_input(
        "Temperature Decay",
        min_value=0.85,
        max_value=0.99,
        value=DEFAULT_PARAMS["temp_decay"],
        step=0.01,
        help="Temperature decay rate per step"
    )
    st.markdown("</div>", unsafe_allow_html=True)
    
    st.markdown("<div class='sidebar-param'>", unsafe_allow_html=True)
    st.session_state.params["num_steps"] = st.number_input(
        "Number of Steps",
        min_value=50,
        max_value=1000,
        value=DEFAULT_PARAMS["num_steps"],
        step=50,
        help="Number of annealing steps"
    )
    st.markdown("</div>", unsafe_allow_html=True)
    
    st.markdown("<div class='sidebar-param'>", unsafe_allow_html=True)
    st.session_state.params["num_repeats"] = st.number_input(
        "Number of Repeats",
        min_value=1,
        max_value=50,
        value=DEFAULT_PARAMS["num_repeats"],
        step=1,
        help="Number of annealing repeats"
    )
    st.markdown("</div>", unsafe_allow_html=True)
    
    # 其他参数
    st.markdown("<div class='sidebar-param'>", unsafe_allow_html=True)
    st.session_state.params["seed"] = st.number_input(
        "Random Seed",
        min_value=1,
        max_value=9999,
        value=DEFAULT_PARAMS["seed"],
        step=1,
        help="Random seed for reproducibility"
    )
    st.markdown("</div>", unsafe_allow_html=True)

# ====================== 主界面 ======================
st.markdown("<div class='main-container'>", unsafe_allow_html=True)

# 标题区域
st.markdown("""
<div class="title-section">
    <h1>pegLIT</h1>
    <p>Automatically identify non-interfering nucleotide linkers between a pegRNA and 3' motif</p>
</div>
""", unsafe_allow_html=True)

# 输入区域
st.markdown("<div class='input-card'>", unsafe_allow_html=True)
st.markdown("<h3 style='font-size:1.5rem; font-weight:600; margin-bottom:16px;'>Sequence Input</h3>", unsafe_allow_html=True)

# 序列输入表格
col_add, col_import, col_export = st.columns([1, 1, 1])
with col_add:
    if st.button("➕ Add New Row", key="add_row", class_="btn-secondary"):
        st.session_state.rows.append(DEFAULT_SEQ.copy())
        st.rerun()

with col_import:
    uploaded_file = st.file_uploader(
        "Import CSV",
        type="csv",
        label_visibility="collapsed",
        key="csv_upload"
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
            
            st.session_state.rows = []
            for _, row in df.iterrows():
                st.session_state.rows.append(row.to_dict())
            st.success("✅ CSV imported successfully!")
            st.rerun()
        except Exception as e:
            st.error(f"❌ CSV import failed: {str(e)}")

with col_export:
    def export_to_csv():
        df = pd.DataFrame(st.session_state.rows)
        buffer = BytesIO()
        df.to_csv(buffer, index=False)
        buffer.seek(0)
        return base64.b64encode(buffer.getvalue()).decode()
    
    if st.session_state.rows:
        csv_base64 = export_to_csv()
        st.download_button(
            label="📤 Export CSV",
            data=base64.b64decode(csv_base64),
            file_name="peglit_sequences.csv",
            mime="text/csv",
            class_="btn-secondary"
        )

# 渲染输入表格
st.markdown("<table class='input-table'>", unsafe_allow_html=True)
st.markdown("""
<thead>
    <tr>
        <th>Spacer</th>
        <th>Scaffold</th>
        <th>Template</th>
        <th>PBS</th>
        <th>Linker Pattern</th>
        <th>Motif</th>
    </tr>
</thead>
<tbody>
""", unsafe_allow_html=True)

for idx, row in enumerate(st.session_state.rows):
    st.markdown("<tr>", unsafe_allow_html=True)
    
    # Spacer列
    st.markdown("<td>", unsafe_allow_html=True)
    st.text_input(
        label=f"spacer_{idx}",
        value=row["spacer"],
        label_visibility="collapsed",
        key=f"spacer_{idx}"
    )
    st.markdown("</td>", unsafe_allow_html=True)
    
    # Scaffold列
    st.markdown("<td>", unsafe_allow_html=True)
    st.text_input(
        label=f"scaffold_{idx}",
        value=row["scaffold"],
        label_visibility="collapsed",
        key=f"scaffold_{idx}"
    )
    st.markdown("</td>", unsafe_allow_html=True)
    
    # Template列
    st.markdown("<td>", unsafe_allow_html=True)
    st.text_input(
        label=f"template_{idx}",
        value=row["template"],
        label_visibility="collapsed",
        key=f"template_{idx}"
    )
    st.markdown("</td>", unsafe_allow_html=True)
    
    # PBS列
    st.markdown("<td>", unsafe_allow_html=True)
    st.text_input(
        label=f"pbs_{idx}",
        value=row["pbs"],
        label_visibility="collapsed",
        key=f"pbs_{idx}"
    )
    st.markdown("</td>", unsafe_allow_html=True)
    
    # Linker列
    st.markdown("<td>", unsafe_allow_html=True)
    st.text_input(
        label=f"linker_{idx}",
        value=row["linker"],
        label_visibility="collapsed",
        disabled=st.session_state.calculated,  # 计算后锁定
        key=f"linker_{idx}"
    )
    st.markdown("</td>", unsafe_allow_html=True)
    
    # Motif列
    st.markdown("<td>", unsafe_allow_html=True)
    st.text_input(
        label=f"motif_{idx}",
        value=row["motif"],
        label_visibility="collapsed",
        key=f"motif_{idx}"
    )
    st.markdown("</td>", unsafe_allow_html=True)
    
    st.markdown("</tr>", unsafe_allow_html=True)

st.markdown("</tbody></table>", unsafe_allow_html=True)

# 计算按钮
st.markdown("<div style='text-align:center; margin:24px 0;'>", unsafe_allow_html=True)
if st.button("START CALCULATION", key="start_calc", class_="btn-primary"):
    with st.spinner("🔄 Running pegLIT calculation... Please wait"):
        st.session_state.results = {}
        st.session_state.calculated = True
        
        try:
            for i, row in enumerate(st.session_state.rows):
                # 获取输入序列（兜底+格式化）
                spacer = row.get("spacer", "").upper().strip()
                scaffold = row.get("scaffold", DEFAULT_SEQ["scaffold"]).upper().strip()
                template = row.get("template", "").upper().strip()
                pbs = row.get("pbs", "").upper().strip()
                motif = row.get("motif", "").upper().strip()
                linker_pattern = row.get("linker", DEFAULT_SEQ["linker"]).upper().strip()
                
                # 空值校验
                if not all([spacer, scaffold, template, pbs, motif]):
                    st.session_state.results[i] = {
                        "status": "error",
                        "message": "Missing required sequence (Spacer/Template/PBS/Motif)"
                    }
                    continue
                
                # 调用核心计算逻辑（复用第一个参考的peglit_min）
                result = peglit_min.pegLIT(
                    seq_spacer=spacer,
                    seq_scaffold=scaffold,
                    seq_template=template,
                    seq_pbs=pbs,
                    seq_motif=motif,
                    linker_pattern=linker_pattern,
                    **st.session_state.params
                )
                
                # 解析结果（兼容官网输出格式）
                row_result = {"status": "success", "data": {}}
                
                # 1. 提取最优linker（对齐官网优先级）
                best_linker = DEFAULT_SEQ["linker"]
                if isinstance(result, dict):
                    best_linker = result.get("best_linker", result.get("linker", best_linker))
                    row_result["data"]["linker"] = best_linker
                    
                    # 2. 提取二级结构（RNAfold预测）
                    full_seq = f"{scaffold}{best_linker}{motif}"
                    (ss, mfe) = RNA.fold(full_seq)
                    row_result["data"]["secondary_structure"] = ss
                    row_result["data"]["mfe"] = mfe
                    
                    # 3. 提取候选列表
                    row_result["data"]["candidates"] = result.get("candidates", [])
                    
                elif isinstance(result, list) and len(result) > 0:
                    best_linker = result[0] if isinstance(result[0], str) else result[0].get("linker", best_linker)
                    row_result["data"]["linker"] = best_linker
                    
                    # 补充二级结构
                    full_seq = f"{scaffold}{best_linker}{motif}"
                    (ss, mfe) = RNA.fold(full_seq)
                    row_result["data"]["secondary_structure"] = ss
                    row_result["data"]["mfe"] = mfe
                    
                elif isinstance(result, str):
                    best_linker = result
                    row_result["data"]["linker"] = best_linker
                    
                    # 补充二级结构
                    full_seq = f"{scaffold}{best_linker}{motif}"
                    (ss, mfe) = RNA.fold(full_seq)
                    row_result["data"]["secondary_structure"] = ss
                    row_result["data"]["mfe"] = mfe
                
                # 更新session state
                st.session_state.rows[i]["linker"] = best_linker
                st.session_state.results[i] = row_result
            
            st.success("✅ Calculation completed successfully!")
            st.rerun()
            
        except Exception as e:
            st.error(f"❌ Calculation failed: {str(e)}")
            st.exception(e)
            st.session_state.calculated = False

st.markdown("</div>", unsafe_allow_html=True)
st.markdown("</div>", unsafe_allow_html=True)  # 关闭input-card

# 结果展示区域
if st.session_state.calculated and st.session_state.results:
    st.markdown("<div class='result-card'>", unsafe_allow_html=True)
    st.markdown("<h3>Calculation Results</h3>", unsafe_allow_html=True)
    
    for idx, result in st.session_state.results.items():
        st.markdown(f"<h4 style='margin:20px 0 10px;'>Row {idx+1}</h4>", unsafe_allow_html=True)
        
        if result["status"] == "error":
            st.error(f"❌ {result['message']}")
            continue
        
        # 核心结果展示
        result_data = result["data"]
        col1, col2 = st.columns([1, 2])
        
        with col1:
            st.markdown("<strong>Optimal Linker:</strong>", unsafe_allow_html=True)
            st.markdown(f"<div style='font-size:1.2rem; font-weight:600; color:#2563eb;'>{result_data['linker']}</div>", unsafe_allow_html=True)
            
            st.markdown("<strong>Minimum Free Energy (MFE):</strong>", unsafe_allow_html=True)
            st.markdown(f"<div style='font-size:1rem;'>{result_data['mfe']:.2f} kcal/mol</div>", unsafe_allow_html=True)
        
        with col2:
            st.markdown("<strong>RNA Secondary Structure:</strong>", unsafe_allow_html=True)
            st.markdown(f"<div class='rna-structure'>{result_data['secondary_structure']}</div>", unsafe_allow_html=True)
        
        # 候选linker列表（如果有）
        if "candidates" in result_data and result_data["candidates"]:
            st.markdown("<strong>Top Candidate Linkers:</strong>", unsafe_allow_html=True)
            candidates_df = pd.DataFrame(result_data["candidates"])
            st.dataframe(
                candidates_df,
                use_container_width=True,
                hide_index=True
            )
        
        st.markdown("<hr style='margin:20px 0; border:1px solid #e5e7eb;'>", unsafe_allow_html=True)
    
    st.markdown("</div>", unsafe_allow_html=True)

st.markdown("</div>", unsafe_allow_html=True)  # 关闭main-container
