import streamlit as st
st.set_page_config(page_title="pegLIT", layout="wide")

import numpy as np
import pandas as pd
import RNA
import peglit_min

if "rows" not in st.session_state:
    st.session_state.rows = [
        {"spacer": "", "scaffold": "GTTTTAG...", "template": "", "pbs": "", "linker": "NNNNNNNN", "motif": "tevopreQ₁"}
    ]

st.markdown("""
<style>
/* 全局 */
* {
    font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif !important;
    box-sizing: border-box;
}
body {
    background-color: #fff;
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
}

/* 表头 */
.table-header {
    display: grid;
    grid-template-columns: 1fr 1.5fr 1.5fr 1fr 1fr 1.5fr;
    gap: 0;
    background: #fff;
    padding: 1rem 0;
    font-weight: 500;
    font-size: 1.1rem;
    border-bottom: 1px solid #e5e7eb;
    background-image: linear-gradient(to right, #e5e7eb 1px, transparent 1px);
    background-size: calc(100% / 6) 100%;
}
.table-header > div {
    padding: 0 1rem;
    text-align: left;
}

/* 输入行 */
.table-input-row {
    display: grid;
    grid-template-columns: 1fr 1.5fr 1.5fr 1fr 1fr 1.5fr;
    gap: 0;
    border-bottom: 1px solid #e5e7eb;
    align-items: center;
    background-image: linear-gradient(to right, #e5e7eb 1px, transparent 1px);
    background-size: calc(100% / 6) 100%;
}
.table-input-row input {
    width: 100%;
    border: none !important;
    outline: none !important;
    font-size: 1rem;
    padding: 0.8rem 1rem;
    background: transparent !important;
}
.table-input-row input:focus {
    background: #f3f4f6 !important;
}

/* 按钮行 */
.action-row {
    display: flex;
    gap: 10px;
    padding: 0.8rem 1rem;
    align-items: center;
}

/* 加号按钮 */
.add-btn {
    width: 34px;
    height: 34px;
    border-radius: 50%;
    border: 1px solid #d1d5db;
    background: white;
    font-size: 16px;
    color: #6b7280;
    display: flex;
    align-items: center;
    justify-content: center;
    cursor: pointer;
}

/* 小上传图标 */
.csv-icon {
    width: 20px;
    height: 20px;
    color: #6b7280;
    cursor: pointer;
}

/* 🔥 彻底隐藏所有原生上传区域 */
.stFileUploader {
    display: none !important;
    visibility: hidden !important;
    height: 0 !important;
    width: 0 !important;
    margin: 0 !important;
    padding: 0 !important;
}

/* START按钮 */
.stButton>button[kind="primary"] {
    background-color: #3b82f6 !important;
    color: white !important;
    border: none !important;
    border-radius: 8px !important;
    padding: 0.8rem 2.5rem !important;
    font-size: 1.4rem !important;
}

/* 隐藏默认元素 */
#MainMenu, footer, header {visibility: hidden;}
</style>
""", unsafe_allow_html=True)

# 标题
st.markdown("<h1>pegLIT</h1>", unsafe_allow_html=True)
st.markdown("""
<div class="subtitle">
Automatically identify non-interfering nucleotide<br>
linkers between a pegRNA and 3' motif.
</div>
""", unsafe_allow_html=True)

# 表格
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

updated_rows = []
for idx, row in enumerate(st.session_state.rows):
    st.markdown("<div class='table-input-row'>", unsafe_allow_html=True)
    cols = st.columns([1, 1.5, 1.5, 1, 1, 1.5])
    updated_row = {
        "spacer": cols[0].text_input("", value=row["spacer"], label_visibility="collapsed", key=f"sp_{idx}"),
        "scaffold": cols[1].text_input("", value=row["scaffold"], label_visibility="collapsed", key=f"sc_{idx}"),
        "template": cols[2].text_input("", value=row["template"], label_visibility="collapsed", key=f"t_{idx}"),
        "pbs": cols[3].text_input("", value=row["pbs"], label_visibility="collapsed", key=f"p_{idx}"),
        "linker": cols[4].text_input("", value=row["linker"], label_visibility="collapsed", key=f"lk_{idx}"),
        "motif": cols[5].text_input("", value=row["motif"], label_visibility="collapsed", key=f"m_{idx}"),
    }
    updated_rows.append(updated_row)
    st.markdown("</div>", unsafe_allow_html=True)

# 加号 + 小图标
st.markdown("<div class='action-row'>", unsafe_allow_html=True)

col1, col2 = st.columns([1, 10])
with col1:
    if st.button("⊕", key="add", use_container_width=False):
        st.session_state.rows.append({"spacer":"","scaffold":"","template":"","pbs":"","linker":"NNNNNNNN","motif":""})
        st.rerun()

with col2:
    # 小箭头图标
    st.markdown("""
    <svg class="csv-icon" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2">
      <path d="M21 15v4a2 2 0 01-2 2H5a2 2 0 01-2-2v-4"></path>
      <path d="M7 10l5 5 5-5"></path>
      <path d="M12 15V3"></path>
    </svg>
    """, unsafe_allow_html=True)

    # 隐藏上传，但功能保留
    file = st.file_uploader("", type="csv", label_visibility="collapsed", key="csv")
    if file:
        df = pd.read_csv(file)
        df.columns = ["spacer","scaffold","template","pbs","linker","motif"]
        st.session_state.rows = df.to_dict("records")
        st.rerun()

st.markdown("</div></div>", unsafe_allow_html=True)

# START
if st.button("START", type="primary"):
    try:
        for i, r in enumerate(updated_rows):
            res = peglit_min.pegLIT(
                seq_spacer=r["spacer"],
                seq_scaffold=r["scaffold"],
                seq_template=r["template"],
                seq_pbs=r["pbs"],
                seq_motif=r["motif"],
                linker_pattern=r["linker"]
            )
            updated_rows[i]["linker"] = res.iloc[0]['linker']
        st.session_state.rows = updated_rows
        st.success("Done")
        st.rerun()
    except Exception as e:
        st.error(e)

st.session_state.rows = updated_rows
