# 必须放在最开头
import streamlit as st
st.set_page_config(page_title="pegLIT", layout="wide")

import pandas as pd
import numpy as np
from peglit_min import pegLIT  # 直接调用你的算法包

# ================== 全局样式（1:1复刻官网） ==================
st.markdown("""
<style>
/* 全局重置 */
* {
    margin: 0;
    padding: 0;
    box-sizing: border-box;
    font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, "Helvetica Neue", Arial, sans-serif;
}
body {
    background-color: #ffffff;
    color: #212529;
}
/* 标题样式 */
h1 {
    text-align: center;
    font-size: 5rem;
    font-weight: 700;
    margin: 2rem 0 1rem;
    letter-spacing: -0.05em;
}
/* 副标题样式 */
.subtitle {
    text-align: center;
    font-size: 1.75rem;
    color: #6c757d;
    line-height: 1.6;
    margin-bottom: 2.5rem;
}
.subtitle a {
    color: #6c757d;
    text-decoration: none;
    border-bottom: 1px solid #6c757d;
}
.subtitle a:hover {
    color: #415E9B;
    border-bottom-color: #415E9B;
}
/* 表格容器 */
.table-wrapper {
    max-width: 1200px;
    margin: 0 auto;
    padding: 0 1rem;
    border: 1px solid #dee2e6;
    border-radius: 8px;
    overflow: hidden;
}
/* 表头 */
.table-header {
    display: grid;
    grid-template-columns: 1fr 1.5fr 1.5fr 1fr 1fr 1.5fr;
    background-color: #f8f9fa;
    padding: 1rem;
    font-weight: 600;
    font-size: 1.1rem;
    border-bottom: 1px solid #dee2e6;
}
/* 表行 */
.table-row {
    display: grid;
    grid-template-columns: 1fr 1.5fr 1.5fr 1fr 1fr 1.5fr;
    padding: 0.75rem 1rem;
    border-bottom: 1px solid #dee2e6;
    align-items: center;
}
.table-row:last-child {
    border-bottom: none;
}
/* 输入框 */
.table-row input {
    width: 100%;
    border: none;
    outline: none;
    font-size: 1rem;
    padding: 0.5rem;
    background-color: transparent;
}
.table-row input:focus {
    background-color: #f8f9fa;
    border-radius: 4px;
}
/* 操作按钮区 */
.action-bar {
    max-width: 1200px;
    margin: 0.5rem auto 2rem;
    padding: 0 1rem;
    display: flex;
    gap: 0.5rem;
}
.action-btn {
    background: none;
    border: 1px solid #dee2e6;
    border-radius: 4px;
    width: 36px;
    height: 36px;
    font-size: 1.25rem;
    cursor: pointer;
    display: flex;
    align-items: center;
    justify-content: center;
}
.action-btn:hover {
    background-color: #f8f9fa;
    border-color: #415E9B;
}
/* START按钮 */
.start-btn-container {
    text-align: center;
    margin: 2rem 0;
}
.stButton>button {
    background-color: #415E9B;
    color: white;
    border: none;
    border-radius: 4px;
    padding: 0.75rem 2rem;
    font-size: 1.25rem;
    font-weight: 600;
    cursor: pointer;
}
.stButton>button:hover {
    background-color: #324b7a;
}
/* 隐藏Streamlit默认元素 */
#MainMenu {visibility: hidden;}
footer {visibility: hidden;}
header {visibility: hidden;}
</style>
""", unsafe_allow_html=True)

# ================== 页面标题（复刻官网） ==================
st.markdown("<h1>pegLIT</h1>", unsafe_allow_html=True)
st.markdown("""
<div class="subtitle">
Automatically identify non-interfering nucleotide<br>
linkers between a pegRNA and 3' motif.
<br><br>
<a href="#">Learn more...</a>
</div>
""", unsafe_allow_html=True)

# ================== 初始化会话状态（多行输入） ==================
if "rows" not in st.session_state:
    st.session_state.rows = [
        {"Spacer": "", "Scaffold": "", "Template": "", "PBS": "", "Linker Pattern": "NNNNNNNN", "Motif": ""}
    ]

# ================== 操作按钮（+ 添加行 / ↑ 导入CSV） ==================
st.markdown("<div class='action-bar'>", unsafe_allow_html=True)
col1, col2 = st.columns([1, 1])
with col1:
    if st.button("➕", key="add_row", help="Add row"):
        st.session_state.rows.append(
            {"Spacer": "", "Scaffold": "", "Template": "", "PBS": "", "Linker Pattern": "NNNNNNNN", "Motif": ""}
        )
        st.rerun()
with col2:
    uploaded_file = st.file_uploader("⬆️", type="csv", key="import_csv", help="Import CSV", label_visibility="collapsed")
    if uploaded_file is not None:
        df = pd.read_csv(uploaded_file)
        st.session_state.rows = df.to_dict("records")
        st.rerun()
st.markdown("</div>", unsafe_allow_html=True)

# ================== 输入表格（复刻官网） ==================
st.markdown("<div class='table-wrapper'>", unsafe_allow_html=True)
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

# 表行
updated_rows = []
for i, row in enumerate(st.session_state.rows):
    cols = st.columns([1, 1.5, 1.5, 1, 1, 1.5])
    updated_row = {}
    updated_row["Spacer"] = cols[0].text_input(f"spacer_{i}", value=row["Spacer"], label_visibility="collapsed")
    updated_row["Scaffold"] = cols[1].text_input(f"scaffold_{i}", value=row["Scaffold"], label_visibility="collapsed")
    updated_row["Template"] = cols[2].text_input(f"template_{i}", value=row["Template"], label_visibility="collapsed")
    updated_row["PBS"] = cols[3].text_input(f"pbs_{i}", value=row["PBS"], label_visibility="collapsed")
    updated_row["Linker Pattern"] = cols[4].text_input(f"linker_{i}", value=row["Linker Pattern"], label_visibility="collapsed")
    updated_row["Motif"] = cols[5].text_input(f"motif_{i}", value=row["Motif"], label_visibility="collapsed")
    updated_rows.append(updated_row)
st.session_state.rows = updated_rows
st.markdown("</div>", unsafe_allow_html=True)

# ================== START按钮（复刻官网） ==================
st.markdown("<div class='start-btn-container'>", unsafe_allow_html=True)
if st.button("START"):
    df_input = pd.DataFrame(st.session_state.rows)
    if df_input.isnull().values.any() or (df_input == "").values.any():
        st.error("❌ Please fill in all fields!")
        st.stop()

    with st.spinner("🔬 Calculating..."):
        try:
            results = []
            for _, row in df_input.iterrows():
                result = pegLIT(
                    seq_spacer=row["Spacer"],
                    seq_scaffold=row["Scaffold"],
                    seq_template=row["Template"],
                    seq_pbs=row["PBS"],
                    seq_motif=row["Motif"],
                    linker_pattern=row["Linker Pattern"]
                )
                results.append(result)

            final_result = pd.concat(results, ignore_index=True)
            st.success("✅ Calculation complete!")
            st.dataframe(final_result, use_container_width=True)

            csv = final_result.to_csv(index=False)
            st.download_button(
                label="📥 Download Results",
                data=csv,
                file_name="peglit_result.csv",
                mime="text/csv"
            )
        except Exception as e:
            st.error(f"❌ Error: {str(e)}")
            st.exception(e)
st.markdown("</div>", unsafe_allow_html=True)
