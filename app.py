import streamlit as st
import pandas as pd

st.set_page_config(page_title="pegLIT", layout="wide")

# 初始化行数据
if "rows" not in st.session_state:
    st.session_state.rows = [
        {"spacer": "", "scaffold": "", "template": "", "pbs": "", "linker": "", "motif": ""}
    ]

# 样式
st.markdown("""
<style>
/* 按钮透明无边框 */
div[data-testid="stButton"] button {
    background: transparent !important;
    border: none !important;
    box-shadow: none !important;
    padding: 0 !important;
    margin: 0 !important;
}
/* 隐藏页面多余元素 */
#MainMenu, footer, header { visibility: hidden; }
</style>
""", unsafe_allow_html=True)

# 标题
st.markdown("<h1 style='text-align:center'>pegLIT</h1>", unsafe_allow_html=True)

# 表格开始
st.markdown("<div class='table-card'>", unsafe_allow_html=True)

# 表头
st.markdown("""
<div style='display:grid; grid-template-columns:1fr 1.5fr 1.5fr 1fr 1fr 1.5fr; padding:1rem 0; border-bottom:1px solid #eee; font-weight:500;'>
    <div style='padding:0 1rem'>Spacer</div>
    <div style='padding:0 1rem'>Scaffold</div>
    <div style='padding:0 1rem'>Template</div>
    <div style='padding:0 1rem'>PBS</div>
    <div style='padding:0 1rem'>Linker Pattern</div>
    <div style='padding:0 1rem'>Motif</div>
</div>
""", unsafe_allow_html=True)

# 输入行
updated_rows = []
for idx, row in enumerate(st.session_state.rows):
    cols = st.columns([1, 1.5, 1.5, 1, 1, 1.5])

    updated_row = {
        "spacer": cols[0].text_input(f"spacer_{idx}", value=row["spacer"], label_visibility="collapsed"),
        "scaffold": cols[1].text_input(f"scaffold_{idx}", value=row["scaffold"], label_visibility="collapsed"),
        "template": cols[2].text_input(f"template_{idx}", value=row["template"], label_visibility="collapsed"),
        "pbs": cols[3].text_input(f"pbs_{idx}", value=row["pbs"], label_visibility="collapsed"),
        
        # 🔥 Linker 只读，不可输入
        "linker": cols[4].text_input(
            f"linker_{idx}", 
            value=row["linker"], 
            label_visibility="collapsed",
            disabled=True
        ),
        
        "motif": cols[5].text_input(f"motif_{idx}", value=row["motif"], label_visibility="collapsed"),
    }
    updated_rows.append(updated_row)

st.session_state.rows = updated_rows

# 按钮行：加号 + 上传图标
btn_col1, btn_col2 = st.columns([1, 1])

with btn_col1:
    # 加号添加行
    if st.button("⊕", key="add_row"):
        st.session_state.rows.append({
            "spacer": "", "scaffold": "", "template": "", "pbs": "", "linker": "", "motif": ""
        })
        st.rerun()

with btn_col2:
    # 上传图标（箭头）
    st.markdown("""
    <svg width="24" height="24" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2">
      <path d="M21 15v4a2 2 0 0 1-2 2H5a2 2 0 0 1-2-2v-4"></path>
      <path d="M7 10l5 5 5-5"></path>
      <path d="M12 15V3"></path>
    </svg>
    """, unsafe_allow_html=True)

    # 点击弹出上传
    if st.button("Select CSV", key="open_upload", label_visibility="collapsed"):
        st.session_state.show_upload = True

# 隐藏式上传
if "show_upload" in st.session_state and st.session_state.show_upload:
    file = st.file_uploader("Upload CSV", type="csv", label_visibility="collapsed")
    if file:
        df = pd.read_csv(file)
        df.columns = ["spacer", "scaffold", "template", "pbs", "linker", "motif"]
        for i, r in df.iterrows():
            if i < len(st.session_state.rows):
                st.session_state.rows[i]["spacer"]   = str(r["spacer"])
                st.session_state.rows[i]["scaffold"] = str(r["scaffold"])
                st.session_state.rows[i]["template"] = str(r["template"])
                st.session_state.rows[i]["pbs"]      = str(r["pbs"])
                st.session_state.rows[i]["linker"]   = str(r["linker"])
                st.session_state.rows[i]["motif"]    = str(r["motif"])
        st.session_state.show_upload = False
        st.rerun()

st.markdown("</div>", unsafe_allow_html=True)

# START 按钮 + 运行中提示
if st.button("START", type="primary"):
    with st.spinner("🔄 Running... Please wait"):
        # 这里放你的计算逻辑
        # 示例：简单赋值
        for i in range(len(st.session_state.rows)):
            st.session_state.rows[i]["linker"] = "RESULT_" + str(i+1)
        st.success("✅ Done")
        st.rerun()
