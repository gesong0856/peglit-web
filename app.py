import streamlit as st
st.set_page_config(page_title="pegLIT", layout="wide")
# 强制清除缓存
st.cache_data.clear()
st.cache_resource.clear()

import pandas as pd
import RNA
import peglit_min
import random
import numpy as np

# ====================== 1. 初始化（仅保留初始1行，对齐官网8个N默认值） ======================
DEFAULT_SEQ = {
    "spacer": "",
    "scaffold": "GTTTCAGAGCTATGCTGGAAACAGCATAGCAAGTTGAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGC",
    "template": "",
    "pbs": "",
    "linker": "NNNNNNNN",  # 严格对齐官网8个N默认值
    "motif": "",
    "unique_seed": 2020  # 每行独立seed，避免linker重复
}

if "rows" not in st.session_state:
    # 初始仅1行，满足需求
    st.session_state.rows = [DEFAULT_SEQ.copy()]

# ====================== 2. 全局样式（完全保留你的原始样式，无任何修改） ======================
st.markdown("""
<style>
/* 全局重置 */
* {
    font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif !important;
    box-sizing: border-box;
    margin: 0;
    padding: 0;
}
body {
    background-color: #ffffff;
}

/* 标题 */
h1 {
    text-align: center;
    font-size: 3rem;
    font-weight: 700;
    margin: 2rem 0 0.5rem !important;
    color: #1f2937;
}
.subtitle {
    text-align: center;
    font-size: 1.1rem;
    color: #1e40af; /* 深蓝 */
    margin-bottom: 2rem;
    line-height: 1.6;
}

/* 表格卡片 */
.table-card {
    max-width: 1200px;
    margin: 0 auto 1rem;
    border: 1px solid #e5e7eb;
    border-radius: 8px;
    overflow: hidden;
    background: white;
    box-shadow: 0 2px 8px rgba(0,0,0,0.1);
}

/* 表头 */
.table-header {
    display: grid;
    grid-template-columns: 1fr 1.5fr 1.5fr 1fr 1fr 1.5fr;
    gap: 0;
    background: #fff;
    padding: 1.25rem 0;
    font-weight: 500;
    font-size: 1.2rem;
    border-bottom: 1px solid #e5e7eb;
}
.table-header > div {
    padding: 0 1rem;
    text-align: left;
    line-height: 1.5;
}

/* 输入行 */
.table-input-row {
    display: grid;
    grid-template-columns: 1fr 1.5fr 1.5fr 1fr 1fr 1.5fr;
    gap: 0;
    border-bottom: 1px solid #e5e7eb;
    align-items: center;
}

/* 输入框样式（绑定初始值+Linker只读） */
.table-input-row input {
    width: 100%;
    border: none !important;
    outline: none !important;
    font-size: 1rem;
    padding: 0.75rem 1rem;
    background: transparent !important;
    line-height: 1.5;
}
.table-input-row input:disabled {
    background: #f3f4f6 !important;
    color: #1e40af !important; /* 深蓝 */
    cursor: not-allowed !important;
    opacity: 1 !important;
}

/* 按钮行（对齐+可点击） */
.action-row {
    display: flex;
    align-items: center;
    gap: 12px;
    padding: 12px 16px;
    height: 48px;
}

/* 加号按钮 */
.add-btn {
    width: 48px;
    height: 48px;
    border-radius: 12px;
    border: 1px solid #d1d5db;
    background: white;
    font-size: 20px;
    color: #374151;
    display: flex;
    align-items: center;
    justify-content: center;
    cursor: pointer;
    transition: all 0.2s;
}
.add-btn:hover {
    border-color: #3b82f6;
    color: #3b82f6;
    background: #f3f4f6;
}

/* 上传按钮（直接用文件框做按钮，100%可点击） */
div[data-testid="stFileUploader"] {
    width: 48px;
    height: 48px;
    border-radius: 12px;
    border: 1px solid #d1d5db;
    background: white;
    overflow: hidden;
    position: relative;
    cursor: pointer;
    transition: all 0.2s;
}
div[data-testid="stFileUploader"]:hover {
    border-color: #3b82f6;
    background: #f3f4f6;
}
div[data-testid="stFileUploader"] > div {
    display: none !important;
}
div[data-testid="stFileUploader"]::after {
    content: "⬆️";
    position: absolute;
    top: 50%;
    left: 50%;
    transform: translate(-50%, -50%);
    font-size: 24px;
    color: #374151;
    pointer-events: none;
}
div[data-testid="stFileUploader"]::before {
    content: "Import CSV";
    position: absolute;
    bottom: 120%;
    left: 50%;
    transform: translateX(-50%);
    background: #1f2937;
    color: white;
    padding: 6px 10px;
    border-radius: 6px;
    font-size: 12px;
    white-space: nowrap;
    opacity: 0;
    visibility: hidden;
    transition: opacity 0.2s ease-in-out;
    z-index: 999;
    pointer-events: none;
}
div[data-testid="stFileUploader"]:hover::before {
    opacity: 1;
    visibility: visible;
}

/* START按钮（蓝色） */
.stButton>button[kind="primary"] {
    background-color: #2563eb;
    color: white !important;
    border: none !important;
    border-radius: 12px !important;
    padding: 0.8rem 2.5rem !important;
    font-size: 1.2rem !important;
}
.stButton>button[kind="primary"]:hover {
    background-color: #1d4ed8;
    color: white !important;
}

/* 隐藏默认元素 */
#MainMenu, footer, header {visibility: hidden;}
</style>
""", unsafe_allow_html=True)

# ====================== 3. 页面标题（保留原始样式） ======================
st.markdown("<h1>pegLIT</h1>", unsafe_allow_html=True)
st.markdown("""
<div class="subtitle">
Automatically identify non-interfering nucleotide linkers between a pegRNA and 3' motif.
</div>
""", unsafe_allow_html=True)

# ====================== 4. 表格渲染（仅初始1行，保留布局） ======================
st.markdown("<div class='table-card'>", unsafe_allow_html=True)

st.markdown("""
<div class="table-header">
    <div>Spacer</div>
    <div>Scaffold</div>
    <div>Template</div>
    <div>PBS</div>
    <div>Linker Pattern</div>
    <div>Motif</div>
</div>
""", unsafe_allow_html=True)

for idx, row in enumerate(st.session_state.rows):
    st.markdown("<div class='table-input-row'>", unsafe_allow_html=True)
    cols = st.columns([1, 1.5, 1.5, 1, 1, 1.5])
    
    cols[0].text_input(
        label=f"spacer_{idx}",
        value=row["spacer"],
        label_visibility="collapsed",
        key=f"spacer_{idx}"
    )
    cols[1].text_input(
        label=f"scaffold_{idx}",
        value=row["scaffold"],
        label_visibility="collapsed",
        key=f"scaffold_{idx}"
    )
    cols[2].text_input(
        label=f"template_{idx}",
        value=row["template"],
        label_visibility="collapsed",
        key=f"template_{idx}"
    )
    cols[3].text_input(
        label=f"pbs_{idx}",
        value=row["pbs"],
        label_visibility="collapsed",
        key=f"pbs_{idx}"
    )
    # Linker只读逻辑：仅计算出有效结果后禁用（保留8个N默认值可编辑）
    linker_disabled = len(str(row["linker"]).strip()) > 0 and row["linker"] != "NNNNNNNN"
    cols[4].text_input(
        label=f"linker_{idx}",
        value=row["linker"],
        label_visibility="collapsed",
        disabled=linker_disabled,
        key=f"linker_{idx}"
    )
    cols[5].text_input(
        label=f"motif_{idx}",
        value=row["motif"],
        label_visibility="collapsed",
        key=f"motif_{idx}"
    )
    st.markdown("</div>", unsafe_allow_html=True)

# ====================== 5. 操作按钮行（支持新增行，保留原始逻辑） ======================
st.markdown("<div class='action-row'>", unsafe_allow_html=True)

if st.button("⊕", key="add_row", help="Add new row"):
    new_row = DEFAULT_SEQ.copy()
    # 新增行生成唯一seed（避免和已有行重复）
    new_row["unique_seed"] = max([r["unique_seed"] for r in st.session_state.rows]) + 100
    st.session_state.rows.append(new_row)
    st.rerun()

# CSV导入（优化序列清洗，对齐官网逻辑）
uploaded_file = st.file_uploader("Upload CSV", type="csv", label_visibility="collapsed", key="csv_upload")
if uploaded_file is not None:
    try:
        df = pd.read_csv(uploaded_file)
        df.columns = df.columns.str.lower()
        required_cols = ["spacer", "scaffold", "template", "pbs", "linker", "motif"]
        for col in required_cols:
            if col not in df.columns:
                df[col] = DEFAULT_SEQ[col]
        df = df[required_cols]
        
        # 重置行数据，导入多少行就生成多少行
        st.session_state.rows = []
        base_seed = 2020
        # 清洗导入序列（大写+去空格，对齐官网）
        for i, (_, row) in enumerate(df.iterrows()):
            row_dict = row.to_dict()
            for k in row_dict:
                row_dict[k] = str(row_dict[k]).upper().strip() if pd.notna(row_dict[k]) else DEFAULT_SEQ[k]
            # 为每行分配唯一seed
            row_dict["unique_seed"] = base_seed + i * 100
            st.session_state.rows.append(row_dict)
        st.success("✅ CSV imported successfully!")
        st.rerun()
    except Exception as e:
        st.error(f"❌ CSV import failed: {str(e)}")

st.markdown("</div></div>", unsafe_allow_html=True)

# ====================== 6. 核心工具函数（官网逻辑+去重机制） ======================
def validate_nucleotide_seq(seq):
    """官网核心：碱基序列合法性校验"""
    valid_chars = set("ATCGN")
    seq = seq.upper().strip()
    if not seq:
        return True, ""
    if set(seq) - valid_chars:
        return False, f"Invalid characters: {set(seq) - valid_chars}. Only A/T/C/G/N are allowed."
    return True, ""

def preprocess_for_peglit(spacer, scaffold, template, pbs, motif, linker_pattern, seed=2020):
    """官网完整预处理逻辑"""
    # 1. 重置随机种子（每行唯一seed，核心去重）
    random.seed(seed)
    np.random.seed(seed)
    
    # 2. 序列标准化（大写+去空格）
    spacer = spacer.upper().strip()
    scaffold = scaffold.upper().strip()
    template = template.upper().strip()
    pbs = pbs.upper().strip()
    motif = motif.upper().strip()
    linker_pattern = linker_pattern.upper().strip() or "NNNNNNNN"  # 兜底8个N
    
    # 3. T/U替换（RNA计算统一用U，结果返回时转回T）
    spacer = spacer.replace("T", "U")
    scaffold = scaffold.replace("T", "U")
    template = template.replace("T", "U")
    pbs = pbs.replace("T", "U")
    motif = motif.replace("T", "U")
    linker_pattern = linker_pattern.replace("T", "U")
    
    # 4. AC阈值修正（比例→绝对数，官网核心）
    ac_thresh_original = 0.5
    ac_thresh = ac_thresh_original * len(linker_pattern)
    
    return {
        "spacer": spacer,
        "scaffold": scaffold,
        "template": template,
        "pbs": pbs,
        "motif": motif,
        "linker_pattern": linker_pattern,
        "ac_thresh": ac_thresh,
        "seed": seed
    }

def get_unique_linker(row_idx, base_linker, used_linkers, max_retries=5):
    """确保linker不重复：如果重复则重新计算"""
    current_linker = base_linker
    retry_count = 0
    # 获取当前行的唯一seed
    current_seed = st.session_state.rows[row_idx]["unique_seed"]
    
    while current_linker in used_linkers and retry_count < max_retries:
        retry_count += 1
        # 重新生成seed，重新计算
        new_seed = current_seed + retry_count * 10
        random.seed(new_seed)
        np.random.seed(new_seed)
        
        # 重新提取当前行序列
        r = st.session_state.rows[row_idx]
        spacer = r.get("spacer", "").upper().strip()
        scaffold = r.get("scaffold", DEFAULT_SEQ["scaffold"]).upper().strip()
        template = r.get("template", "").upper().strip()
        pbs = r.get("pbs", "").upper().strip()
        motif = r.get("motif", "").upper().strip()
        linker_pattern = r.get("linker", DEFAULT_SEQ["linker"]).upper().strip()
        
        # 重新预处理
        preprocessed = preprocess_for_peglit(
            spacer=spacer,
            scaffold=scaffold,
            template=template,
            pbs=pbs,
            motif=motif,
            linker_pattern=linker_pattern,
            seed=new_seed
        )
        
        # 重新调用算法
        try:
            result = peglit_min.pegLIT(
                seq_spacer=preprocessed["spacer"].replace("U", "T"),
                seq_scaffold=preprocessed["scaffold"].replace("U", "T"),
                seq_template=preprocessed["template"].replace("U", "T"),
                seq_pbs=preprocessed["pbs"].replace("U", "T"),
                seq_motif=preprocessed["motif"].replace("U", "T"),
                linker_pattern=preprocessed["linker_pattern"].replace("U", "T"),
                ac_thresh=preprocessed["ac_thresh"],
                u_thresh=3,
                n_thresh=3,
                topn=100,
                epsilon=1e-2,
                num_repeats=10,
                num_steps=250,
                temp_init=0.15,
                temp_decay=0.95,
                bottleneck=1,
                seed=new_seed,
                sequences_to_avoid=None
            )
            # 解析新结果
            if isinstance(result, str):
                current_linker = result.strip().replace("U", "T")
            elif isinstance(result, list) and len(result) > 0:
                if isinstance(result[0], dict):
                    current_linker = result[0].get("linker", "").strip().replace("U", "T")
                else:
                    current_linker = str(result[0]).strip().replace("U", "T")
            elif isinstance(result, dict):
                current_linker = result.get("linker", "").strip().replace("U", "T")
        except Exception as e:
            st.warning(f"Row {row_idx+1}: Retry {retry_count} failed - {str(e)}")
            current_linker = f"NNNNNNNN_{row_idx}_{retry_count}"  # 兜底唯一值
    
    return current_linker

# ====================== 7. START按钮（核心：多行列独立计算+linker去重） ======================
st.markdown("<div class='start-btn-container'>", unsafe_allow_html=True)
if st.button("START", type="primary"):
    with st.spinner("🔄 Running... Please wait"):
        try:
            total_rows = len(st.session_state.rows)
            used_linkers = set()  # 记录已使用的linker，避免重复
            success_count = 0
            
            for i, r in enumerate(st.session_state.rows):
                st.write(f"📌 Processing Row {i+1}/{total_rows}")
                
                # 1. 提取原始输入
                spacer = r.get("spacer", "")
                scaffold = r.get("scaffold", DEFAULT_SEQ["scaffold"])
                template = r.get("template", "")
                pbs = r.get("pbs", "")
                motif = r.get("motif", "")
                linker_pattern = r.get("linker", DEFAULT_SEQ["linker"])
                unique_seed = r["unique_seed"]  # 每行唯一seed
                
                # 2. 序列合法性校验
                all_seq = spacer + scaffold + template + pbs + motif + linker_pattern
                is_valid, err_msg = validate_nucleotide_seq(all_seq)
                if not is_valid:
                    st.error(f"Row {i+1}: ❌ Invalid sequence - {err_msg}")
                    st.session_state.rows[i]["linker"] = f"NNNNNNNN_{i}"  # 唯一兜底值
                    used_linkers.add(st.session_state.rows[i]["linker"])
                    continue
                
                # 3. 核心序列非空校验
                core_seqs = [spacer, scaffold, template, pbs, motif]
                if not all(core_seqs):
                    st.warning(f"Row {i+1}: ⚠️ Core sequences cannot be empty!")
                    st.session_state.rows[i]["linker"] = f"NNNNNNNN_{i}"  # 唯一兜底值
                    used_linkers.add(st.session_state.rows[i]["linker"])
                    continue
                
                # 4. 官网预处理（使用每行唯一seed）
                preprocessed = preprocess_for_peglit(
                    spacer=spacer,
                    scaffold=scaffold,
                    template=template,
                    pbs=pbs,
                    motif=motif,
                    linker_pattern=linker_pattern,
                    seed=unique_seed
                )
                
                # 5. 调用算法（捕获解包错误）
                base_linker = ""
                try:
                    result = peglit_min.pegLIT(
                        seq_spacer=preprocessed["spacer"].replace("U", "T"),
                        seq_scaffold=preprocessed["scaffold"].replace("U", "T"),
                        seq_template=preprocessed["template"].replace("U", "T"),
                        seq_pbs=preprocessed["pbs"].replace("U", "T"),
                        seq_motif=preprocessed["motif"].replace("U", "T"),
                        linker_pattern=preprocessed["linker_pattern"].replace("U", "T"),
                        ac_thresh=preprocessed["ac_thresh"],
                        u_thresh=3,
                        n_thresh=3,
                        topn=100,
                        epsilon=1e-2,
                        num_repeats=10,
                        num_steps=250,
                        temp_init=0.15,
                        temp_decay=0.95,
                        bottleneck=1,
                        seed=preprocessed["seed"],
                        sequences_to_avoid=None
                    )
                    # 解析基础结果
                    if isinstance(result, str):
                        base_linker = result.strip().replace("U", "T")
                    elif isinstance(result, list) and len(result) > 0:
                        if isinstance(result[0], dict):
                            base_linker = result[0].get("linker", "").strip().replace("U", "T")
                        else:
                            base_linker = str(result[0]).strip().replace("U", "T")
                    elif isinstance(result, dict):
                        base_linker = result.get("linker", "").strip().replace("U", "T")
                except ValueError as e:
                    if "not enough values to unpack" in str(e):
                        st.error(f"Row {i+1}: ❌ No valid linker found (filter too strict)")
                        base_linker = f"NNNNNNNN_{i}"  # 唯一兜底值
                    else:
                        raise e
                except Exception as e:
                    st.error(f"Row {i+1}: ❌ Calculation error - {str(e)}")
                    base_linker = f"NNNNNNNN_{i}"  # 唯一兜底值
                
                # 6. 核心：确保linker唯一（去重逻辑）
                final_linker = get_unique_linker(i, base_linker, used_linkers)
                
                # 7. 更新状态并记录已使用的linker
                st.session_state.rows[i]["linker"] = final_linker
                used_linkers.add(final_linker)
                
                # 8. 结果提示
                if final_linker != f"NNNNNNNN_{i}" and not final_linker.startswith("NNNNNNNN_"):
                    st.success(f"Row {i+1}: ✅ Unique linker generated - {final_linker}")
                    success_count += 1
                else:
                    st.warning(f"Row {i+1}: ⚠️ Unique default linker - {final_linker}")
            
            # 最终汇总提示
            st.success(f"🎉 All calculations completed! Successfully generated {success_count}/{total_rows} unique linkers.")
            st.rerun()
        except Exception as e:
            st.error(f"❌ Global error: {str(e)}")
            st.exception(e)

st.markdown("</div>", unsafe_allow_html=True)
