import streamlit as st
st.set_page_config(page_title="pegLIT", layout="wide")

import numpy as np
import pandas as pd
import RNA
import peglit_min

# ====================== 1. 初始化会话状态（你的默认序列） ======================
if "rows" not in st.session_state:
    st.session_state.rows = [
        {
            "spacer": "ACCCTGCCTTGCTAAGGCCA",
            "scaffold": "GTTTCAGAGCTATGCTGGAAACAGCATAGCAAGTTGAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGC",
            "template": "TTTCAATGTTCTCTCTGATcGCGTATTTGCATGCcctCTGGCAGATTTCTGTgATgTCAGCACCACTGAAACCTTGcGTaTATTTGGCAAGAGCATTCAGaTCTACgTCCTTGGCCACAGGTGACTTtctGAGGCAAGCTTTGAAGATCTGCAGgcgAGATTGATCATCAGGCAGtGGgATGTAGATAAGCTGATCAAGACGaCCTGG",
            "pbs": "CCTTAGCAAG",
            "linker": "NNNNNNNN",
            "motif": "TTGACGCGGTTCTATCTAGTTACGCGTTAAACCAACTAGAAA"
        }
    ]

if "show_upload" not in st.session_state:
    st.session_state.show_upload = False

# ====================== 2. 全局样式 ======================
st.markdown("""
<style>
/* 全局 */
* {
    font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif !important;
    box-sizing: border-box;
}
body {
    background-color: #ffffff;
}

/* 标题 */
h1 {
    text-align: center;
    font-size: 2.5rem;
    margin: 2rem 0 0.5rem !important;
    color: #1f2937;
}
.subtitle {
    text-align: center;
    font-size: 1rem;
    color: #6b7280;
    margin-bottom: 2rem;
}

/* 表格 */
.table-card {
    max-width: 1400px;
    margin: 0 auto 1rem;
    border: 1px solid #e5e7eb;
    border-radius: 8px;
    background: white;
}
.table-header {
    display: grid;
    grid-template-columns: 1fr 1.8fr 2.5fr 1fr 1fr 1.8fr;
    padding: 1rem 0;
    font-weight: 500;
    border-bottom: 1px solid #eee;
}
.table-header > div {
    padding: 0 1rem;
}
.table-input-row {
    display: grid;
    grid-template-columns: 1fr 1.8fr 2.5fr 1fr 1fr 1.8fr;
    border-bottom: 1px solid #eee;
}

/* 输入框 */
.table-input-row input {
    border: none !important;
    outline: none !important;
    padding: 0.75rem 1rem;
    background: transparent !important;
}
.table-input-row input:disabled {
    background: #f3f4f6 !important;
    color: #444 !important;
}

/* 按钮行 */
.action-row {
    display: flex;
    align-items: center;
    gap: 12px;
    padding: 12px 16px;
}

/* 上传按钮样式 */
.upload-button-wrapper {
    display: flex;
    align-items: center;
    justify-content: center;
    width: 44px;
    height: 44px;
    border-radius: 8px;
    cursor: pointer;
}
.upload-button-wrapper:hover {
    background: #f3f4f6;
}

/* 隐藏原生上传 */
div[data-testid="stFileUploader"] {
    display: none !important;
}
#MainMenu, footer, header {visibility: hidden;}
</style>
""", unsafe_allow_html=True)

# ====================== 标题 ======================
st.markdown("<h1>pegLIT</h1>", unsafe_allow_html=True)
st.markdown("""
<div class="subtitle">
Automatically identify non-interfering nucleotide linkers between a pegRNA and 3' motif.
</div>
""", unsafe_allow_html=True)

# ====================== 表格 ======================
st.markdown("<div class='table-card'>", unsafe_allow_html=True)

st.markdown("""
<div class='table-header'>
    <div>Spacer</div>
    <div>Scaffold</div>
    <div>Template</div>
    <div>PBS</div>
    <div>Linker Pattern</div>
    <div>Motif</div>
</div>
""", unsafe_allow_html=True)

updated_rows = []
for idx, row in enumerate(st.session_state.rows):
    st.markdown("<div class='table-input-row'>", unsafe_allow_html=True)
    cols = st.columns([1, 1.8, 2.5, 1, 1, 1.8])

    updated_row = {
        "spacer": cols[0].text_input(f"spacer_{idx}", value=row["spacer"], label_visibility="collapsed"),
        "scaffold": cols[1].text_input(f"scaffold_{idx}", value=row["scaffold"], label_visibility="collapsed"),
        "template": cols[2].text_input(f"template_{idx}", value=row["template"], label_visibility="collapsed"),
        "pbs": cols[3].text_input(f"pbs_{idx}", value=row["pbs"], label_visibility="collapsed"),
        "linker": cols[4].text_input(f"linker_{idx}", value=row["linker"], label_visibility="collapsed", disabled=True),
        "motif": cols[5].text_input(f"motif_{idx}", value=row["motif"], label_visibility="collapsed"),
    }
    updated_rows.append(updated_row)
    st.markdown("</div>", unsafe_allow_html=True)

# ====================== 按钮：加号 + 上传 ======================
st.markdown("<div class='action-row'>", unsafe_allow_html=True)

# 加号
if st.button("⊕", key="add", help="Add row"):
    st.session_state.rows.append({
        "spacer": "", "scaffold": "", "template": "", "pbs": "", "linker": "NNNNNNNN", "motif": ""
    })
    st.rerun()

# 上传按钮（完美对齐、有提示、可点击）
st.markdown("""
<div class="upload-button-wrapper" title="Import CSV">
    <svg width="20" height="20" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2">
        <path d="M21 15v4a2 2 0 0 1-2 2H5a2 2 0 0 1-2-2v-4"></path>
        <path d="M7 10l5-5 5 5"></path>
        <path d="M12 15V3"></path>
    </svg>
</div>
""", unsafe_allow_html=True)

# 点击触发上传
if st.button("", key="upload_trigger", help="Import CSV"):
    st.session_state.show_upload = True
    st.rerun()

st.markdown("</div></div>", unsafe_allow_html=True)

# ====================== 上传功能 ======================
if st.session_state.show_upload:
    f = st.file_uploader("CSV", type="csv")
    if f:
        df = pd.read_csv(f)
        df.columns = ["spacer", "scaffold", "template", "pbs", "linker", "motif"]
        for i, r in df.iterrows():
            if i < len(st.session_state.rows):
                st.session_state.rows[i]["spacer"] = r["spacer"]
                st.session_state.rows[i]["scaffold"] = r["scaffold"]
                st.session_state.rows[i]["template"] = r["template"]
                st.session_state.rows[i]["pbs"] = r["pbs"]
                st.session_state.rows[i]["linker"] = r["linker"]
                st.session_state.rows[i]["motif"] = r["motif"]
        st.session_state.show_upload = False
        st.rerun()

# ====================== START ======================
if st.button("START", type="primary", use_container_width=False):
    with st.spinner("🔄 Running..."):
        try:
            for i, r in enumerate(updated_rows):
                out = peglit_min.pegLIT(
                    seq_spacer=r["spacer"],
                    seq_scaffold=r["scaffold"],
                    seq_template=r["template"],
                    seq_pbs=r["pbs"],
                    seq_motif=r["motif"],
                    linker_pattern=r["linker"]
                )
                updated_rows[i]["linker"] = out.iloc[0]['linker']
            st.session_state.rows = updated_rows
            st.success("✅ Done")
            st.rerun()
        except Exception as e:
            st.error(f"Error: {e}")

st.session_state.rows = updated_rows
