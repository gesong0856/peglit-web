import streamlit as st
st.set_page_config(page_title="pegLIT", layout="wide")
# 强制清除缓存
st.cache_data.clear()
st.cache_resource.clear()

import pandas as pd
import RNA
import peglit_min

# ====================== 1. 初始化（增强健壮性） ======================
DEFAULT_SEQ = {
    "spacer": "ACCCTGCCTTGCTAAGGCCA",
    "scaffold": "GTTTCAGAGCTATGCTGGAAACAGCATAGCAAGTTGAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGC",
    "template": "TTTCAATGTTCTCTCTGATcGCGTATTTGCATGCcctCTGGCAGATTTCTGTgATgTCAGCACCACTGAAACCTTGcGTaTATTTGGCAAGAGCATTCAGaTCTACgTCCTTGGCCACAGGTGACTTtctGAGGCAAGCTTTGAAGATCTGCAGgcgAGATTGATCATCAGGCAGtGGgATGTAGATAAGCTGATCAAGACGaCCTGG",
    "pbs": "CCTTAGCAAG",
    "linker": "NNNNNNNN",
    "motif": "TTGACGCGGTTCTATCTAGTTACGCGTTAAACCAACTAGAAA"
}

if "rows" not in st.session_state:
    st.session_state.rows = [DEFAULT_SEQ.copy()]
# 新增：初始化存储计算结果的状态，用于额外展示
if "calculated_linkers" not in st.session_state:
    st.session_state.calculated_linkers = []

# ====================== 2. 全局样式 ======================
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

# ====================== 3. 页面标题 ======================
st.markdown("<h1>pegLIT</h1>", unsafe_allow_html=True)
st.markdown("""
<div class="subtitle">
Automatically identify non-interfering nucleotide linkers between a pegRNA and 3' motif.
</div>
""", unsafe_allow_html=True)

# ====================== 4. 表格渲染 ======================
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
    cols[4].text_input(
        label=f"linker_{idx}",
        value=row["linker"],
        label_visibility="collapsed",
        disabled=True,
        key=f"linker_{idx}"
    )
    cols[5].text_input(
        label=f"motif_{idx}",
        value=row["motif"],
        label_visibility="collapsed",
        key=f"motif_{idx}"
    )
    st.markdown("</div>", unsafe_allow_html=True)

# ====================== 5. 操作按钮行 ======================
st.markdown("<div class='action-row'>", unsafe_allow_html=True)

if st.button("⊕", key="add_row", help="Add new row"):
    st.session_state.rows.append(DEFAULT_SEQ.copy())
    st.rerun()

# 上传按钮 + 修复CSV导入逻辑
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
            if i < len(st.session_state.rows):
                st.session_state.rows[i] = row_dict
            else:
                st.session_state.rows.append(row_dict)
        st.success("✅ CSV imported successfully!")
        st.rerun()
    except Exception as e:
        st.error(f"❌ CSV import failed: {str(e)}")

st.markdown("</div></div>", unsafe_allow_html=True)

# ====================== 6. START按钮（核心修复：兼容所有返回格式） ======================
st.markdown("<div class='start-btn-container'>", unsafe_allow_html=True)
if st.button("START", type="primary"):
    with st.spinner("🔄 Running... Please wait"):
        try:
            # 新增：清空历史计算结果
            st.session_state.calculated_linkers = []
            for i, r in enumerate(st.session_state.rows):
                # 关键修复：每个键都做兜底，避免缺失
                spacer = r.get("spacer", DEFAULT_SEQ["spacer"])
                scaffold = r.get("scaffold", DEFAULT_SEQ["scaffold"])
                template = r.get("template", DEFAULT_SEQ["template"])
                pbs = r.get("pbs", DEFAULT_SEQ["pbs"])  # 重点：pbs键兜底
                motif = r.get("motif", DEFAULT_SEQ["motif"])
                linker = r.get("linker", DEFAULT_SEQ["linker"])
                
                # 调用工具
                result = peglit_min.pegLIT(
                    seq_spacer=spacer,
                    seq_scaffold=scaffold,
                    seq_template=template,
                    seq_pbs=pbs,
                    seq_motif=motif,
                    linker_pattern=linker,
                    ac_thresh=0.5,
                    u_thresh=3,
                    n_thresh=3,
                    topn=100,
                    epsilon=1e-2,
                    num_repeats=10,
                    num_steps=250,
                    temp_init=0.15,
                    temp_decay=0.9,
                    bottleneck=1,
                    seed=2020,
                    sequences_to_avoid=None
                )
                
                # ========== 核心修复：兼容字符串/列表/字典所有返回格式 ==========
                new_linker = DEFAULT_SEQ["linker"]  # 默认值
                if isinstance(result, str):
                    # 情况1：算法直接返回linker字符串（如"ATCGATCG"）
                    new_linker = result
                elif isinstance(result, list) and len(result) > 0:
                    # 情况2：返回列表 → 兼容列表内是字典/字符串
                    if isinstance(result[0], dict):
                        new_linker = result[0].get("linker", DEFAULT_SEQ["linker"])
                    else:
                        # 列表内是字符串（如["ATCGATCG"]）
                        new_linker = result[0]
                elif isinstance(result, dict):
                    # 情况3：返回字典
                    new_linker = result.get("linker", DEFAULT_SEQ["linker"])
                # 情况4：其他格式（None/空）→ 保留默认值
                
                # 更新linker并提示
                st.session_state.rows[i]["linker"] = new_linker
                # 新增：将计算结果存入状态，用于额外展示
                st.session_state.calculated_linkers.append({
                    "row_num": i+1,
                    "linker": new_linker
                })
                if new_linker == DEFAULT_SEQ["linker"]:
                    st.warning(f"Row {i+1}: No valid linker result, keeping default (NNNNNNNN).")
                else:
                    st.success(f"Row {i+1}: Linker updated to {new_linker}")
            
            st.success("✅ Calculation completed!")
            st.rerun()  # 强制刷新界面显示新结果
        except Exception as e:
            st.error(f"❌ Error: {str(e)}")
            st.exception(e)

st.markdown("</div>", unsafe_allow_html=True)
