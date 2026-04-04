# 必须放在最开头
import streamlit as st
st.set_page_config(page_title="pegLIT", layout="wide")

import pandas as pd
import numpy as np
import RNA
import peglit_min  # 你的核心算法文件

# ================== 全局样式（复刻官网） ==================
st.markdown("""
<style>
/* 全局字体和背景 */
body {
    font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, "Helvetica Neue", Arial, sans-serif;
    background-color: #ffffff;
}
/* 标题样式 */
h1 {
    text-align: center;
    font-size: 4rem;
    font-weight: 700;
    color: #212529;
    margin-top: 2rem;
    margin-bottom: 1rem;
}
/* 副标题样式 */
.subtitle {
    text-align: center;
    font-size: 1.75rem;
    color: #6c757d;
    margin-bottom: 2rem;
    line-height: 1.5;
}
/* 表格容器 */
.table-container {
    max-width: 1200px;
    margin: 0 auto;
    padding: 0 1rem;
}
/* 按钮样式 */
.stButton>button {
    background-color: #415E9B;
    color: white;
    border: none;
    border-radius: 4px;
    padding: 0.75rem 2rem;
    font-size: 1.25rem;
    font-weight: 600;
    margin: 2rem auto;
    display: block;
    cursor: pointer;
}
.stButton>button:hover {
    background-color: #324b7a;
}
/* 操作按钮（+和↑） */
.action-btn {
    background: none;
    border: none;
    font-size: 1.5rem;
    cursor: pointer;
    padding: 0.25rem;
    margin: 0 0.25rem;
}
</style>
""", unsafe_allow_html=True)

# ================== 页面标题（复刻官网） ==================
st.markdown("<h1>pegLIT</h1>", unsafe_allow_html=True)
st.markdown("""
<div class="subtitle">
Automatically identify non-interfering nucleotide<br>
linkers between a pegRNA and 3' motif.
<br><br>
<a href="#" style="color: #6c757d; text-decoration: none;">Learn more...</a>
</div>
""", unsafe_allow_html=True)

# ================== 初始化会话状态（多行输入） ==================
if "rows" not in st.session_state:
    # 初始1行数据
    st.session_state.rows = [{
        "Spacer": "",
        "Scaffold": "",
        "Template": "",
        "PBS": "",
        "Linker Pattern": "NNNNNNNN",
        "Motif": ""
    }]

# ================== 操作按钮区（添加行/导入CSV） ==================
col1, col2 = st.columns([1, 1])
with col1:
    if st.button("➕", key="add_row", help="添加一行"):
        # 新增空行
        st.session_state.rows.append({
            "Spacer": "",
            "Scaffold": "",
            "Template": "",
            "PBS": "",
            "Linker Pattern": "NNNNNNNN",
            "Motif": ""
        })
        st.rerun()
with col2:
    uploaded_file = st.file_uploader("⬆️", type="csv", key="import_csv", help="导入CSV文件", label_visibility="collapsed")
    if uploaded_file is not None:
        # 读取CSV并更新行数据
        df = pd.read_csv(uploaded_file)
        st.session_state.rows = df.to_dict("records")
        st.rerun()

# ================== 输入表格（复刻官网） ==================
st.markdown("<div class='table-container'>", unsafe_allow_html=True)

# 表头
cols = st.columns(6)
headers = ["Spacer", "Scaffold", "Template", "PBS", "Linker Pattern", "Motif"]
for col, header in zip(cols, headers):
    col.write(f"**{header}**")

# 输入行
updated_rows = []
for i, row in enumerate(st.session_state.rows):
    cols = st.columns(6)
    updated_row = {}
    updated_row["Spacer"] = cols[0].text_input(f"spacer_{i}", value=row["Spacer"], label_visibility="collapsed")
    updated_row["Scaffold"] = cols[1].text_input(f"scaffold_{i}", value=row["Scaffold"], label_visibility="collapsed")
    updated_row["Template"] = cols[2].text_input(f"template_{i}", value=row["Template"], label_visibility="collapsed")
    updated_row["PBS"] = cols[3].text_input(f"pbs_{i}", value=row["PBS"], label_visibility="collapsed")
    updated_row["Linker Pattern"] = cols[4].text_input(f"linker_{i}", value=row["Linker Pattern"], label_visibility="collapsed")
    updated_row["Motif"] = cols[5].text_input(f"motif_{i}", value=row["Motif"], label_visibility="collapsed")
    updated_rows.append(updated_row)

# 更新会话状态
st.session_state.rows = updated_rows
st.markdown("</div>", unsafe_allow_html=True)

# ================== 运行按钮（复刻官网） ==================
if st.button("START"):
    # 校验输入
    df_input = pd.DataFrame(st.session_state.rows)
    if df_input.isnull().values.any() or (df_input == "").values.any():
        st.error("❌ 请填写所有输入项！")
        st.stop()

    with st.spinner("🔬 正在计算..."):
        try:
            # 批量计算（支持多行）
            results = []
            for _, row in df_input.iterrows():
                result = peglit_min.pegLIT(
                    seq_spacer=row["Spacer"],
                    seq_scaffold=row["Scaffold"],
                    seq_template=row["Template"],
                    seq_pbs=row["PBS"],
                    seq_motif=row["Motif"],
                    linker_pattern=row["Linker Pattern"],
                    ac_thresh=0.5,
                    u_thresh=3,
                    n_thresh=3,
                    topn=100
                )
                results.append(result)

            # 合并结果
            final_result = pd.concat(results, ignore_index=True)
            
            # 展示结果
            st.success("✅ 计算完成！")
            st.dataframe(final_result, use_container_width=True)

            # 下载结果
            csv = final_result.to_csv(index=False)
            st.download_button(
                label="📥 下载结果",
                data=csv,
                file_name="peglit_result.csv",
                mime="text/csv"
            )

        except Exception as e:
            st.error(f"❌ 计算出错：{str(e)}")
            st.exception(e)
