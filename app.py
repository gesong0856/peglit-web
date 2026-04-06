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
import string  # 官网依赖：碱基校验

# ====================== 1. 初始化（保留原有逻辑，仅优化空行生成） ======================
DEFAULT_SEQ = {
    "spacer": "",
    "scaffold": "GTTTCAGAGCTATGCTGGAAACAGCATAGCAAGTTGAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGC",
    "template": "",
    "pbs": "",
    "linker": "",  # 取消默认NNNNNNNN，对齐官网空初始值
    "motif": ""
}

if "rows" not in st.session_state:
    st.session_state.rows = [DEFAULT_SEQ.copy()]

# ====================== 2. 全局样式（完全保留你的原始样式，未做任何修改） ======================
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

# ====================== 3. 页面标题（保留你的原始样式） ======================
st.markdown("<h1>pegLIT</h1>", unsafe_allow_html=True)
st.markdown("""
<div class="subtitle">
Automatically identify non-interfering nucleotide linkers between a pegRNA and 3' motif.
</div>
""", unsafe_allow_html=True)

# ====================== 4. 表格渲染（保留原有布局，仅优化linker只读逻辑） ======================
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
    # 优化linker只读逻辑：仅当计算出有效结果后禁用
    linker_disabled = len(str(row["linker"]).strip()) > 0 and row["linker"] != "NNNNNNNN"
    cols[4].text_input(
        label=f"linker_{idx}",
        value=row["linker"],
        label_visibility="collapsed",
        disabled=linker_disabled,  # 计算后只读，未计算可编辑
        key=f"linker_{idx}"
    )
    cols[5].text_input(
        label=f"motif_{idx}",
        value=row["motif"],
        label_visibility="collapsed",
        key=f"motif_{idx}"
    )
    st.markdown("</div>", unsafe_allow_html=True)

# ====================== 5. 操作按钮行（保留原有逻辑） ======================
st.markdown("<div class='action-row'>", unsafe_allow_html=True)

if st.button("⊕", key="add_row", help="Add new row"):
    st.session_state.rows.append(DEFAULT_SEQ.copy())
    st.rerun()

# 上传按钮 + 修复CSV导入逻辑（保留原有逻辑，仅优化空值处理）
uploaded_file = st.file_uploader("Upload CSV", type="csv", label_visibility="collapsed", key="csv_upload")
if uploaded_file is not None:
    try:
        df = pd.read_csv(uploaded_file)
        # 强制对齐列名（忽略大小写+补全缺失列）
        df.columns = df.columns.str.lower()  # 统一列名为小写
        required_cols = ["spacer", "scaffold", "template", "pbs", "linker", "motif"]
        # 补全缺失的列，用默认值填充
        for col in required_cols:
            if col not in df.columns:
                df[col] = DEFAULT_SEQ[col]
        # 只保留需要的列，避免多余列干扰
        df = df[required_cols]
        
        # 导入数据：不足的行覆盖，多余的行新增
        for i, (_, row) in enumerate(df.iterrows()):
            row_dict = row.to_dict()
            # 清洗导入的序列（对齐官网逻辑）
            for k in row_dict:
                row_dict[k] = str(row_dict[k]).upper().strip() if pd.notna(row_dict[k]) else ""
            if i < len(st.session_state.rows):
                st.session_state.rows[i] = row_dict
            else:
                st.session_state.rows.append(row_dict)
        st.success("✅ CSV imported successfully!")
        st.rerun()
    except Exception as e:
        st.error(f"❌ CSV import failed: {str(e)}")

st.markdown("</div></div>", unsafe_allow_html=True)

# ====================== 6. 官网核心工具函数（补全，不改动样式） ======================
def validate_nucleotide_seq(seq):
    """官网核心：碱基序列校验"""
    valid_chars = set("ATCGN")
    seq = seq.upper().strip()
    if not seq:
        return True, ""  # 空序列允许
    if set(seq) - valid_chars:
        return False, f"Invalid characters: {set(seq) - valid_chars}. Only A/T/C/G/N are allowed."
    return True, ""

def preprocess_sequences(spacer, scaffold, template, pbs, motif, linker_pattern, seed=2020):
    """官网完整预处理逻辑"""
    # 1. 随机种子初始化（官网执行顺序）
    random.seed(seed)
    np.random.seed(seed)
    
    # 2. 序列大写+去空格
    spacer = spacer.upper().strip()
    scaffold = scaffold.upper().strip()
    template = template.upper().strip()
    pbs = pbs.upper().strip()
    motif = motif.upper().strip()
    linker_pattern = linker_pattern.upper().strip() or "NNNNNNNN"
    
    # 3. 拼接上下文序列（官网核心：用于二级结构计算）
    seq_pre = spacer + scaffold + template + pbs
    seq_post = motif
    
    # 4. AC阈值修正（比例→绝对数）
    ac_thresh_original = 0.5
    ac_thresh = ac_thresh_original * len(linker_pattern)
    
    # 5. 禁用序列预处理
    sequences_to_avoid = None
    if sequences_to_avoid is not None:
        sequences_to_avoid = set(m.upper() for m in sequences_to_avoid)
    
    # 6. T/U替换（官网隐性处理：RNA计算统一用U）
    seq_pre = seq_pre.replace("T", "U")
    seq_post = seq_post.replace("T", "U")
    linker_pattern = linker_pattern.replace("T", "U")
    
    return {
        "spacer": spacer,
        "scaffold": scaffold,
        "template": template,
        "pbs": pbs,
        "motif": motif,
        "linker_pattern": linker_pattern,
        "seq_pre": seq_pre,
        "seq_post": seq_post,
        "ac_thresh": ac_thresh,
        "sequences_to_avoid": sequences_to_avoid,
        "seed": seed
    }

# ====================== 7. START按钮（补全官网所有核心逻辑，保留样式） ======================
st.markdown("<div class='start-btn-container'>", unsafe_allow_html=True)
if st.button("START", type="primary"):
    with st.spinner("🔄 Running... Please wait"):
        try:
            for i, r in enumerate(st.session_state.rows):
                st.write(f"正在计算 Row {i+1}...")
                
                # 1. 提取原始输入
                spacer = r.get("spacer", "")
                scaffold = r.get("scaffold", DEFAULT_SEQ["scaffold"])
                template = r.get("template", "")
                pbs = r.get("pbs", "")
                motif = r.get("motif", "")
                linker_pattern = r.get("linker", "")
                
                # 2. 官网核心：序列合法性校验
                is_valid, err_msg = validate_nucleotide_seq(spacer + scaffold + template + pbs + motif + linker_pattern)
                if not is_valid:
                    st.error(f"Row {i+1} 序列校验失败：{err_msg}")
                    continue
                
                # 3. 官网完整预处理
                preprocessed = preprocess_sequences(
                    spacer=spacer,
                    scaffold=scaffold,
                    template=template,
                    pbs=pbs,
                    motif=motif,
                    linker_pattern=linker_pattern,
                    seed=2020
                )
                
                # 4. 调用官网算法（参数100%对齐）
                result = peglit_min.pegLIT(
                    seq_spacer=preprocessed["spacer"],
                    seq_scaffold=preprocessed["scaffold"],
                    seq_template=preprocessed["template"],
                    seq_pbs=preprocessed["pbs"],
                    seq_motif=preprocessed["motif"],
                    linker_pattern=preprocessed["linker_pattern"],
                    ac_thresh=preprocessed["ac_thresh"],  # 修正后的AC阈值
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
                    sequences_to_avoid=preprocessed["sequences_to_avoid"]
                )
                
                # 5. 官网结果解析逻辑（兼容所有返回格式）
                new_linker = ""
                if isinstance(result, str):
                    new_linker = result.strip().replace("U", "T")  # RNA→DNA，U转回T
                elif isinstance(result, list) and len(result) > 0:
                    if isinstance(result[0], dict):
                        new_linker = result[0].get("linker", "").strip().replace("U", "T")
                    else:
                        new_linker = str(result[0]).strip().replace("U", "T")
                elif isinstance(result, dict):
                    new_linker = result.get("linker", "").strip().replace("U", "T")
                
                # 6. 兜底处理（无结果时用默认值）
                if not new_linker:
                    new_linker = "NNNNNNNN"
                    st.warning(f"Row {i+1}: No valid linker result, keeping default (NNNNNNNN).")
                else:
                    st.success(f"Row {i+1}: Linker updated to {new_linker}")
                
                # 7. 更新session state（保留原有逻辑）
                st.session_state.rows[i]["linker"] = new_linker
            
            st.success("✅ Calculation completed!")
            st.rerun()
        except Exception as e:
            st.error(f"❌ Error: {str(e)}")
            st.exception(e)

st.markdown("</div>", unsafe_allow_html=True)
