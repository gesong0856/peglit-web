# 【必须放在最最开头！】第一个 Streamlit 命令
import streamlit as st
st.set_page_config(page_title="pegLIT", layout="wide")

# 后续导入
import numpy as np
import pandas as pd
import RNA
import peglit_min

# 新版ViennaRNA默认支持长序列，删除旧版API调用
# RNA.set_parameter("max_length", 2000)

# ================== 会话状态：支持多行输入 ==================
if "rows" not in st.session_state:
    st.session_state.rows = [
        {"spacer": "", "scaffold": "GTTTTAG...", "template": "", "pbs": "", "linker": "NNNNNNNN", "motif": "tevopreQ₁\nCGCGGT..."}
    ]

# ================== 样式：隐藏原生拖拽区，自定义下载图标 ==================
st.markdown("""
<style>
/* 全局重置 */
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
    font-size: 3rem;
    font-weight: 700;
    margin: 2rem 0 0.5rem !important;
    color: #1f2937;
}

/* 副标题 */
.subtitle {
    text-align: center;
    font-size: 1.1rem;
    color: #6b7280;
    margin-bottom: 2rem;
    line-height: 1.6;
}

/* 表格容器 */
.table-card {
    max-width: 1200px;
    margin: 0 auto 1rem;
    border: 1px solid #e5e7eb;
    border-radius: 8px;
    overflow: hidden;
    background: white;
    box-shadow: 0 2px 8px rgba(0,0,0,0.1);
}

/* 表头行：完美居中+单元格内对齐 */
.table-header {
    display: grid;
    grid-template-columns: 1fr 1.5fr 1.5fr 1fr 1fr 1.5fr;
    gap: 0;
    background-color: #ffffff;
    padding: 1.25rem 0;
    font-weight: 500;
    font-size: 1.2rem;
    border-bottom: 1px solid #e5e7eb;
    background-image: linear-gradient(to right, #e5e7eb 1px, transparent 1px);
    background-size: calc(100% / 6) 100%;
    background-repeat: repeat-x;
}
.table-header > div {
    padding: 0 1rem;
    text-align: left;
    line-height: 1.5;
}

/* 输入行：与表头完全同步 */
.table-input-row {
    display: grid;
    grid-template-columns: 1fr 1.5fr 1.5fr 1fr 1fr 1.5fr;
    gap: 0;
    border-bottom: 1px solid #e5e7eb;
    align-items: center;
    background-image: linear-gradient(to right, #e5e7eb 1px, transparent 1px);
    background-size: calc(100% / 6) 100%;
    background-repeat: repeat-x;
}

/* 输入框样式：单元格内对齐 */
.table-input-row input {
    width: 100%;
    border: none !important;
    outline: none !important;
    font-size: 1.1rem;
    padding: 0.75rem 1rem;
    background-color: transparent !important;
    line-height: 1.5;
    text-align: left;
    -webkit-appearance: none;
    appearance: none;
}
.table-input-row input:focus {
    background-color: #f3f4f6 !important;
}

/* 操作按钮行：加号+下载图标 */
.action-row {
    display: grid;
    grid-template-columns: 0.5fr 0.5fr 5fr;
    gap: 0.5rem;
    padding: 0.8rem 1rem;
    align-items: center;
}

/* 圆圈加号按钮 */
.circle-btn {
    width: 36px;
    height: 36px;
    border-radius: 50%;
    border: 1px solid #d1d5db;
    background: white;
    font-size: 18px;
    color: #6b7280;
    cursor: pointer;
    display: flex;
    align-items: center;
    justify-content: center;
}
.circle-btn:hover {
    border-color: #3b82f6;
    color: #3b82f6;
}

/* 下载图标样式 */
.download-icon {
    width: 32px;
    height: 32px;
    fill: #6b7280;
    cursor: pointer;
}
.download-icon:hover {
    fill: #3b82f6;
}

/* 【关键】完全隐藏原生文件上传组件 */
.stFileUploader {
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
    font-size: 1.5rem !important;
    font-weight: 600 !important;
    box-shadow: 0 2px 8px rgba(0,0,0,0.1) !important;
    margin: 0 auto !important;
    display: block !important;
}
.stButton>button[kind="primary"]:hover {
    background-color: #2563eb !important;
}

/* 隐藏Streamlit默认元素 */
#MainMenu, footer, header {visibility: hidden;}
</style>
""", unsafe_allow_html=True)

# ================== 页面标题 ==================
st.markdown("<h1>pegLIT</h1>", unsafe_allow_html=True)
st.markdown("""
<div class="subtitle">
Automatically identify non-interfering nucleotide<br>
linkers between a pegRNA and 3' motif.
</div>
""", unsafe_allow_html=True)

# ================== 表格 ==================
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

# 输入行
updated_rows = []
for idx, row in enumerate(st.session_state.rows):
    st.markdown("<div class='table-input-row'>", unsafe_allow_html=True)
    cols = st.columns([1, 1.5, 1.5, 1, 1, 1.5])
    
    updated_row = {}
    updated_row["spacer"] = cols[0].text_input("", value=row["spacer"], label_visibility="collapsed", key=f"sp_{idx}")
    updated_row["scaffold"] = cols[1].text_input("", value=row["scaffold"], label_visibility="collapsed", key=f"sc_{idx}")
    updated_row["template"] = cols[2].text_input("", value=row["template"], label_visibility="collapsed", key=f"t_{idx}")
    updated_row["pbs"] = cols[3].text_input("", value=row["pbs"], label_visibility="collapsed", key=f"p_{idx}")
    updated_row["linker"] = cols[4].text_input("", value=row["linker"], label_visibility="collapsed", key=f"lk_{idx}")
    updated_row["motif"] = cols[5].text_input("", value=row["motif"], label_visibility="collapsed", key=f"m_{idx}")
    
    updated_rows.append(updated_row)
    st.markdown("</div>", unsafe_allow_html=True)

# 操作按钮行：加号 + 下载图标（隐藏原生上传）
st.markdown("<div class='action-row'>", unsafe_allow_html=True)
col_add, col_csv, _ = st.columns([0.5, 0.5, 5])

# 圆圈加号：添加行
with col_add:
    if st.button("⊕", key="add_row", help="Add row"):
        st.session_state.rows.append({
            "spacer": "", "scaffold": "", "template": "", "pbs": "", "linker": "NNNNNNNN", "motif": ""
        })
        st.rerun()

# 下载图标：绑定文件上传功能，隐藏原生组件
with col_csv:
    # 嵌入下载样式SVG图标
    download_svg = """
    <svg class="download-icon" viewBox="0 0 24 24" xmlns="http://www.w3.org/2000/svg">
        <path d="M21 15V19C21 19.5304 20.7893 20.0391 20.4142 20.4142C20.0391 20.7893 19.5304 21 19 21H5C4.46957 21 3.96086 20.7893 3.58579 20.4142C3.21071 20.0391 3 19.5304 3 19V15M17 8L12 13L7 8M12 13V3" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round"/>
    </svg>
    """
    st.markdown(download_svg, unsafe_allow_html=True)
    # 隐藏式文件上传，点击图标触发
    uploaded_file = st.file_uploader("", type="csv", label_visibility="collapsed", key="csv_upload")
    if uploaded_file is not None:
        df = pd.read_csv(uploaded_file)
        df.columns = ["spacer", "scaffold", "template", "pbs", "linker", "motif"]
        st.session_state.rows = df.to_dict("records")
        st.rerun()

st.markdown("</div></div>", unsafe_allow_html=True)

# ================== START按钮 ==================
st.markdown("<div class='start-btn-container'>", unsafe_allow_html=True)
if st.button("START", type="primary"):
    # 输入校验
    df_input = pd.DataFrame(updated_rows)
    if df_input.isnull().values.any() or (df_input == "").values.any():
        st.error("❌ 请填写所有序列！")
        st.stop()

    with st.spinner("🔬 计算中，请稍候..."):
        try:
            # 批量计算，结果写入Linker Pattern单元格
            for i, row in enumerate(updated_rows):
                result = peglit_min.pegLIT(
                    seq_spacer=row["spacer"],
                    seq_scaffold=row["scaffold"],
                    seq_template=row["template"],
                    seq_pbs=row["pbs"],
                    seq_motif=row["motif"],
                    linker_pattern=row["linker"],
                    ac_thresh=0.5,
                    u_thresh=3,
                    n_thresh=3,
                    topn=1,
                    epsilon=1e-2,
                    num_repeats=10,
                    num_steps=250,
                    temp_init=0.15,
                    temp_decay=0.9,
                    bottleneck=1,
                    seed=2026,
                    sequences_to_avoid=None
                )
                # 最优结果写入单元格
                best_linker = result.iloc[0]['linker']
                updated_rows[i]["linker"] = best_linker

            # 保存结果，刷新界面
            st.session_state.rows = updated_rows
            st.success("✅ 计算完成！结果已显示在Linker Pattern列")
            st.rerun()

        except Exception as e:
            st.error(f"❌ 计算出错：{str(e)}")
            st.exception(e)

st.markdown("</div>", unsafe_allow_html=True)

# 更新会话状态
st.session_state.rows = updated_rows
