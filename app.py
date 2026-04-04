import streamlit as st
st.set_page_config(page_title="pegLIT", layout="wide")

import numpy as np
import pandas as pd
import RNA
import peglit_min

# ====================== 1. 初始化会话状态 ======================
if "rows" not in st.session_state:
    st.session_state.rows = [
        {"spacer": "", "scaffold": "GTTTTAG...", "template": "", "pbs": "", "linker": "NNNNNNNN", "motif": "tevopreQ₁"}
    ]
# 上传弹窗状态控制
if "show_upload_dialog" not in st.session_state:
    st.session_state.show_upload_dialog = False

# ====================== 2. 全局样式（彻底隐藏上传区） ======================
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
    background-image: linear-gradient(to right, #e5e7eb 1px, transparent 1px);
    background-size: calc(100% / 6) 100%;
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
    background-image: linear-gradient(to right, #e5e7eb 1px, transparent 1px);
    background-size: calc(100% / 6) 100%;
    background-repeat: repeat-x;
}

/* 输入框样式（修复空标签+Linker只读样式） */
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
/* Linker 只读样式：灰色背景、不可点击、光标禁用 */
.table-input-row input:disabled {
    color: #6b7280 !important;
    background-color: #f3f4f6 !important;
    cursor: not-allowed !important;
}

/* 操作按钮行 */
.action-row {
    display: flex;
    gap: 12px;
    padding: 0.8rem 1rem;
    align-items: center;
}

/* 加号按钮 */
.add-btn {
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
    transition: all 0.2s;
}
.add-btn:hover {
    border-color: #3b82f6;
    color: #3b82f6;
    background: #f3f4f6;
}

/* 下载图标按钮 */
.download-btn {
    width: 36px;
    height: 36px;
    background: transparent;
    border: none;
    cursor: pointer;
    display: flex;
    align-items: center;
    justify-content: center;
    transition: all 0.2s;
}
.download-btn:hover {
    transform: translateY(-2px);
}

/* 🔥 终极隐藏原生上传区（多层级高优先级，彻底隐藏） */
div[data-testid="stFileUploader"] {
    display: none !important;
    visibility: hidden !important;
    height: 0 !important;
    min-height: 0 !important;
    width: 0 !important;
    min-width: 0 !important;
    margin: 0 !important;
    padding: 0 !important;
    position: absolute !important;
    z-index: -9999 !important;
}
div[data-testid="stFileUploaderDropzone"] {
    display: none !important;
    visibility: hidden !important;
    height: 0 !important;
    min-height: 0 !important;
}
[class*="stFileUploader"] {
    display: none !important;
    visibility: hidden !important;
    height: 0 !important;
    min-height: 0 !important;
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

# ====================== 3. 页面标题 ======================
st.markdown("<h1>pegLIT</h1>", unsafe_allow_html=True)
st.markdown("""
<div class="subtitle">
Automatically identify non-interfering nucleotide<br>
linkers between a pegRNA and 3' motif.
</div>
""", unsafe_allow_html=True)

# ====================== 4. 表格渲染 ======================
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

# 输入行（修复空标签+Linker只读）
updated_rows = []
for idx, row in enumerate(st.session_state.rows):
    st.markdown("<div class='table-input-row'>", unsafe_allow_html=True)
    cols = st.columns([1, 1.5, 1.5, 1, 1, 1.5])
    
    # 给所有输入框加唯一隐形label，彻底消除空标签警告
    updated_row = {
        "spacer": cols[0].text_input(f"spacer_{idx}", value=row["spacer"], label_visibility="collapsed"),
        "scaffold": cols[1].text_input(f"scaffold_{idx}", value=row["scaffold"], label_visibility="collapsed"),
        "template": cols[2].text_input(f"template_{idx}", value=row["template"], label_visibility="collapsed"),
        "pbs": cols[3].text_input(f"pbs_{idx}", value=row["pbs"], label_visibility="collapsed"),
        # Linker Pattern 只读，仅输出，不可手动输入
        "linker": cols[4].text_input(f"linker_{idx}", value=row["linker"], label_visibility="collapsed", disabled=True),
        "motif": cols[5].text_input(f"motif_{idx}", value=row["motif"], label_visibility="collapsed")
    }
    updated_rows.append(updated_row)
    st.markdown("</div>", unsafe_allow_html=True)

# ====================== 5. 操作按钮行（加号+下载图标） ======================
st.markdown("<div class='action-row'>", unsafe_allow_html=True)

# 1. 加号按钮：添加新行
if st.button("⊕", key="add_row_btn", help="Add new row"):
    st.session_state.rows.append({
        "spacer": "", "scaffold": "", "template": "", "pbs": "", "linker": "NNNNNNNN", "motif": ""
    })
    st.rerun()

# 2. 下载图标按钮：点击触发上传（用状态控制，无非法参数）
# 用st.button包裹SVG，仅传合法参数
upload_svg = """
<svg width="24" height="24" viewBox="0 0 24 24" fill="none" xmlns="http://www.w3.org/2000/svg">
  <path d="M21 15v4a2 2 0 0 1-2 2H5a2 2 0 0 1-2-2v-4"></path>
  <path d="M7 10l5 5 5-5"></path>
  <path d="M12 15V3"></path>
</svg>
"""
if st.button(upload_svg, key="download_upload_btn", unsafe_allow_html=True):
    st.session_state.show_upload_dialog = True
    st.rerun()

st.markdown("</div></div>", unsafe_allow_html=True)

# ====================== 6. 隐藏式上传弹窗 ======================
if st.session_state.show_upload_dialog:
    # 仅用合法参数调用st.file_uploader，无label_visibility错误
    uploaded_file = st.file_uploader("Upload CSV", type="csv", key="csv_upload")
    if uploaded_file is not None:
        df = pd.read_csv(uploaded_file)
        # 强制列名匹配，避免导入错误
        df.columns = ["spacer", "scaffold", "template", "pbs", "linker", "motif"]
        # 上传数据自动填入对应行
        for i, row in df.iterrows():
            if i < len(st.session_state.rows):
                st.session_state.rows[i]["spacer"] = str(row["spacer"])
                st.session_state.rows[i]["scaffold"] = str(row["scaffold"])
                st.session_state.rows[i]["template"] = str(row["template"])
                st.session_state.rows[i]["pbs"] = str(row["pbs"])
                st.session_state.rows[i]["linker"] = str(row["linker"])
                st.session_state.rows[i]["motif"] = str(row["motif"])
        # 关闭弹窗，刷新界面
        st.session_state.show_upload_dialog = False
        st.rerun()

# ====================== 7. START按钮+运行中提示 ======================
st.markdown("<div class='start-btn-container'>", unsafe_allow_html=True)
if st.button("START", type="primary"):
    with st.spinner("🔄 Running... Please wait"):
        try:
            # 批量计算pegLIT结果，写入Linker Pattern
            for i, r in enumerate(updated_rows):
                result = peglit_min.pegLIT(
                    seq_spacer=r["spacer"],
                    seq_scaffold=r["scaffold"],
                    seq_template=r["template"],
                    seq_pbs=r["pbs"],
                    seq_motif=r["motif"],
                    linker_pattern=r["linker"]
                )
                # 提取最优结果，写入Linker
                updated_rows[i]["linker"] = result.iloc[0]['linker']
            
            # 更新会话状态，刷新界面
            st.session_state.rows = updated_rows
            st.success("✅ Calculation completed!")
            st.rerun()
        except Exception as e:
            st.error(f"❌ Error: {str(e)}")

st.markdown("</div>", unsafe_allow_html=True)

# ====================== 8. 同步会话状态 ======================
st.session_state.rows = updated_rows
