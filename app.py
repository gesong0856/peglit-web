# 必须放在最开头
import streamlit as st
st.set_page_config(page_title="pegLIT", layout="wide")

import pandas as pd
import numpy as np
from peglit_min import pegLIT

# ==============================================
# 【仅修改这里：官网原版样式 + 按钮位置】
# ==============================================
st.markdown("""
<style>
/* 全局 */
* {
    font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif !important;
}
body {
    background-color: #ffffff;
}

/* 标题 */
h1 {
    text-align: center;
    font-size: 5rem;
    font-weight: 700;
    margin: 2rem 0 0.5rem !important;
    color: #111827;
}

/* 副标题 */
.subtitle {
    text-align: center;
    font-size: 1.3rem;
    color: #6b7280;
    margin-bottom: 3rem;
    line-height: 1.6;
}

/* 表格容器（官网卡片） */
.table-container {
    max-width: 1000px;
    margin: 0 auto 2rem;
    border: 1px solid #e5e7eb;
    border-radius: 10px;
    overflow: hidden;
    background: white;
}

/* 表头 */
.table-header {
    background-color: #f9fafb;
    border-bottom: 1px solid #e5e7eb;
    padding: 14px 20px;
    font-weight: 500;
    font-size: 1rem;
    color: #374151;
}

/* 按钮：圆圈加号 + 上传箭头（官网同款位置） */
.action-row {
    padding: 12px 20px;
    display: flex;
    gap: 12px;
    align-items: center;
    border-bottom: 1px solid #e5e7eb;
}

/* 圆圈加号按钮 */
.circle-btn {
    width: 28px;
    height: 28px;
    border-radius: 50%;
    border: 1px solid #d1d5db;
    background: white;
    font-size: 18px;
    color: #4b5563;
    cursor: pointer;
    display: flex;
    align-items: center;
    justify-content: center;
}
.circle-btn:hover {
    border-color: #2563eb;
    color: #2563eb;
}

/* 上传图标 */
.upload-icon {
    font-size: 20px;
    color: #4b5563;
    cursor: pointer;
}
.upload-icon:hover {
    color: #2563eb;
}

/* 输入行 */
.input-row {
    padding: 10px 20px;
    border-bottom: 1px solid #f3f4f6;
}

/* 隐藏 Streamlit 自带元素 */
#MainMenu, footer, header {visibility: hidden;}

/* START 按钮 */
.start-button {
    display: block;
    margin: 40px auto;
    background-color: #3b82f6;
    color: white;
    border: none;
    border-radius: 8px;
    padding: 14px 36px;
    font-size: 18px;
    font-weight: 500;
    cursor: pointer;
}
.start-button:hover {
    background-color: #2563eb;
}
</style>
""", unsafe_allow_html=True)

# ==============================================
# 标题（官网一模一样）
# ==============================================
st.markdown("<h1>pegLIT</h1>", unsafe_allow_html=True)
st.markdown("""
<div class="subtitle">
Automatically identify non-interfering nucleotide<br>
linkers between a pegRNA and 3' motif.
</div>
""", unsafe_allow_html=True)

# ==============================================
# 初始化多行
# ==============================================
if "rows" not in st.session_state:
    st.session_state.rows = [
        {"Spacer": "", "Scaffold": "", "Template": "", "PBS": "", "Linker Pattern": "NNNNNNNN", "Motif": ""}
    ]

# ==============================================
# 【官网同款界面：表格 + 圆圈加号 + 上传箭头】
# ==============================================
st.markdown("<div class='table-container'>", unsafe_allow_html=True)

# 表头
st.markdown("""
<div class="table-header">
    Spacer  Scaffold  Template  PBS  Linker Pattern  Motif
</div>
""", unsafe_allow_html=True)

# 操作按钮行：圆圈加号 + 上传箭头（官网位置）
st.markdown("<div class='action-row'>", unsafe_allow_html=True)
col1, col2, col3 = st.columns([0.1, 0.1, 1])
with col1:
    if st.button("⊕", key="add", help="Add row"):
        st.session_state.rows.append({
            "Spacer": "", "Scaffold": "", "Template": "", "PBS": "", "Linker Pattern": "NNNNNNNN", "Motif": ""
        })
        st.rerun()
with col2:
    file = st.file_uploader("", type="csv", label_visibility="collapsed")
    if file:
        df = pd.read_csv(file)
        st.session_state.rows = df.to_dict("records")
        st.rerun()
st.markdown("</div>", unsafe_allow_html=True)

# 输入区域（官网同款布局）
new_rows = []
for i, row in enumerate(st.session_state.rows):
    c1, c2, c3, c4, c5, c6 = st.columns([1, 1.5, 1.5, 1, 1, 1.5])
    new_rows.append({
        "Spacer": c1.text_input("", row["Spacer"], key=f"s{i}", label_visibility="collapsed"),
        "Scaffold": c2.text_input("", row["Scaffold"], key=f"sc{i}", label_visibility="collapsed"),
        "Template": c3.text_input("", row["Template"], key=f"t{i}", label_visibility="collapsed"),
        "PBS": c4.text_input("", row["PBS"], key=f"p{i}", label_visibility="collapsed"),
        "Linker Pattern": c5.text_input("", row["Linker Pattern"], key=f"lp{i}", label_visibility="collapsed"),
        "Motif": c6.text_input("", row["Motif"], key=f"m{i}", label_visibility="collapsed"),
    })
st.session_state.rows = new_rows

st.markdown("</div>", unsafe_allow_html=True)

# ==============================================
# START 按钮（官网同款）
# ==============================================
if st.button("START", use_container_width=False, type="primary"):
    inputs = pd.DataFrame(st.session_state.rows)
    if inputs.isna().any().any() or (inputs == "").any().any():
        st.error("Please fill all fields.")
        st.stop()

    with st.spinner("Running..."):
        results = []
        for _, r in inputs.iterrows():
            res = pegLIT(
                seq_spacer=r["Spacer"],
                seq_scaffold=r["Scaffold"],
                seq_template=r["Template"],
                seq_pbs=r["PBS"],
                seq_motif=r["Motif"],
                linker_pattern=r["Linker Pattern"]
            )
            results.append(res)

        final = pd.concat(results, ignore_index=True)
        st.success("Done!")
        st.dataframe(final, use_container_width=True)

        st.download_button("Download CSV", final.to_csv(index=False), "peglit_results.csv")
