import streamlit as st
st.set_page_config(page_title="pegLIT", layout="wide")
# 强制清除缓存
st.cache_data.clear()
st.cache_resource.clear()

import pandas as pd
import RNA
import peglit_min

# ====================== 1. 初始化（只执行一次，不被覆盖） ======================
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
if "show_upload" not in st.session_state:
    st.session_state.show_upload = False

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
    color: #6b7280;
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
    color: #6b7280 !important;
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

/* 上传按钮（直接用图标做按钮，100%可点击） */
.upload-btn {
    width: 48px;
    height: 48px;
    border-radius: 12px;
    border: 1px solid transparent;
    background: #f3f4f6;
    display: flex;
    align-items: center;
    justify-content: center;
    cursor: pointer;
    transition: all 0.2s;
    position: relative;
}
.upload-btn:hover {
    background: #e5e7eb;
}
/* hover提示 */
.upload-btn::after {
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
}
.upload-btn:hover::after {
    opacity: 1;
    visibility: visible;
}

/* 上传图标SVG（内联在按钮中） */
.upload-icon {
    width: 24px;
    height: 24px;
    fill: #374151;
}

/* 隐藏原生上传区 */
div[data-testid="stFileUploader"] {
    display: none !important;
}

/* START按钮 */
.start-btn-container {
    text-align: center;
    margin: 1rem 0 2rem;
}
.stButton>button[kind="primary"] {
    background-color: #3b82f6 !important;
    color: white !important;
    border: none !important;
    border-radius: 8px !important;
    padding: 0.8rem 2.5rem !important;
    font-size: 1.2rem !important;
}
.stButton>button[kind="primary"]:hover {
    background-color: #2563eb !important;
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

# ====================== 4. 表格渲染（正确绑定初始值） ======================
st.markdown("<div class='table-card'>", unsafe_allow_html=True)

# 表头
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

# 输入行（遍历session_state，初始值100%显示）
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

# ====================== 5. 操作按钮行（核心修复：按钮直接带图标，100%可点击） ======================
st.markdown("<div class='action-row'>", unsafe_allow_html=True)

# 1. 加号按钮
if st.button("⊕", key="add_row", help="Add new row"):
    st.session_state.rows.append(DEFAULT_SEQ.copy())
    st.rerun()

# 2. 上传按钮（直接用图标做按钮，无层级冲突，点击秒触发）
upload_icon_svg = """
<svg width="24" height="24" viewBox="0 0 24 24" fill="none" xmlns="http://www.w3.org/2000/svg">
    <path d="M21 15V19C21 19.5304 20.7893 20.0391 20.4142 20.4142C20.0391 20.7893 19.5304 21 19 21H5C4.46957 21 3.96086 20.7893 3.58579 20.4142C3.21071 20.0391 3 19.5304 3 19V15M17 8L12 13L7 8M12 13V3" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round"/>
</svg>
"""
# 按钮直接渲染图标，无层级冲突
if st.button(upload_icon_svg, key="upload_btn", unsafe_allow_html=True):
    st.session_state.show_upload = True
    st.rerun()

st.markdown("</div></div>", unsafe_allow_html=True)

# ====================== 6. 隐藏式上传 ======================
if st.session_state.show_upload:
    uploaded_file = st.file_uploader("Upload CSV", type="csv", label_visibility="collapsed", key="csv_upload")
    if uploaded_file is not None:
        df = pd.read_csv(uploaded_file)
        df.columns = ["spacer", "scaffold", "template", "pbs", "linker", "motif"]
        for i, row in df.iterrows():
            if i < len(st.session_state.rows):
                st.session_state.rows[i] = row.to_dict()
        st.session_state.show_upload = False
        st.rerun()

# ====================== 7. START按钮 ======================
st.markdown("<div class='start-btn-container'>", unsafe_allow_html=True)
if st.button("START", type="primary"):
    with st.spinner("🔄 Running... Please wait"):
        try:
            for i, r in enumerate(st.session_state.rows):
                result = peglit_min.pegLIT(
                    seq_spacer=r["spacer"],
                    seq_scaffold=r["scaffold"],
                    seq_template=r["template"],
                    seq_pbs=r["pbs"],
                    seq_motif=r["motif"],
                    linker_pattern=r["linker"]
                )
                st.session_state.rows[i]["linker"] = result[0]["linker"]
            st.success("✅ Calculation completed!")
            st.rerun()
        except Exception as e:
            st.error(f"❌ Error: {str(e)}")

st.markdown("</div>", unsafe_allow_html=True)
