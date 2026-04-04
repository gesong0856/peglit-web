import streamlit as st
st.set_page_config(page_title="pegLIT", layout="wide")

import numpy as np
import pandas as pd
import RNA
import peglit_min

# ====================== 1. 【核心修复】初始化会话状态（放在最顶部，优先执行） ======================
DEFAULT_ROW = {
    "spacer": "ACCCTGCCTTGCTAAGGCCA",
    "scaffold": "GTTTCAGAGCTATGCTGGAAACAGCATAGCAAGTTGAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGC",
    "template": "TTTCAATGTTCTCTCTGATcGCGTATTTGCATGCcctCTGGCAGATTTCTGTgATgTCAGCACCACTGAAACCTTGcGTaTATTTGGCAAGAGCATTCAGaTCTACgTCCTTGGCCACAGGTGACTTtctGAGGCAAGCTTTGAAGATCTGCAGgcgAGATTGATCATCAGGCAGtGGgATGTAGATAAGCTGATCAAGACGaCCTGG",
    "pbs": "CCTTAGCAAG",
    "linker": "NNNNNNNN",
    "motif": "TTGACGCGGTTCTATCTAGTTACGCGTTAAACCAACTAGAAA"
}

# 初始化会话状态，确保页面加载时就有默认值
if "rows" not in st.session_state:
    st.session_state.rows = [DEFAULT_ROW.copy()]  # 深度拷贝，避免引用问题
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

/* 输入框样式（Linker只读+空标签修复） */
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
.table-input-row input:disabled {
    color: #6b7280 !important;
    background-color: #f3f4f6 !important;
    cursor: not-allowed !important;
}

/* ============== 按钮行100%对齐+hover提示 ============== */
.action-row {
    display: flex;
    align-items: center;
    gap: 16px;
    padding: 12px 16px;
    height: 56px; /* 统一高度，彻底解决错位 */
}

/* 加号按钮（圆角矩形，和图标完全对齐） */
.add-btn {
    width: 56px;
    height: 56px;
    border-radius: 16px;
    border: 1px solid #d1d5db;
    background: white;
    font-size: 24px;
    color: #374151;
    display: inline-flex;
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

/* 下载图标按钮（和加号完全对齐，加hover提示） */
.download-btn {
    width: 56px;
    height: 56px;
    border-radius: 16px;
    border: 1px solid transparent;
    background: transparent;
    display: inline-flex;
    align-items: center;
    justify-content: center;
    cursor: pointer;
    transition: all 0.2s;
    position: relative;
}
.download-btn:hover {
    background: #f3f4f6;
}
/* hover提示框（悬浮显示Import CSV） */
.download-btn::after {
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
.download-btn:hover::after {
    opacity: 1;
    visibility: visible;
}

/* 下载图标SVG样式（还原你要的样式） */
.download-icon {
    width: 32px;
    height: 32px;
    fill: #374151;
}

/* ============== 彻底隐藏原生上传区 ============== */
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

/* 隐藏默认元素 */
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

# ====================== 4. 表格渲染（读取会话状态，初始值生效） ======================
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

# 输入行（遍历会话状态，初始值自动填入）
updated_rows = []
for idx, row in enumerate(st.session_state.rows):
    st.markdown("<div class='table-input-row'>", unsafe_allow_html=True)
    cols = st.columns([1, 1.5, 1.5, 1, 1, 1.5])
    
    # 【关键】从会话状态读取值，填入输入框
    updated_row = {
        "spacer": cols[0].text_input(f"spacer_{idx}", value=row["spacer"], label_visibility="collapsed"),
        "scaffold": cols[1].text_input(f"scaffold_{idx}", value=row["scaffold"], label_visibility="collapsed"),
        "template": cols[2].text_input(f"template_{idx}", value=row["template"], label_visibility="collapsed"),
        "pbs": cols[3].text_input(f"pbs_{idx}", value=row["pbs"], label_visibility="collapsed"),
        # Linker Pattern 只读，仅输出
        "linker": cols[4].text_input(f"linker_{idx}", value=row["linker"], label_visibility="collapsed", disabled=True),
        "motif": cols[5].text_input(f"motif_{idx}", value=row["motif"], label_visibility="collapsed")
    }
    updated_rows.append(updated_row)
    st.markdown("</div>", unsafe_allow_html=True)

# ====================== 5. 操作按钮行（对齐+hover提示+图标还原） ======================
st.markdown("<div class='action-row'>", unsafe_allow_html=True)

# 1. 加号按钮
if st.button("⊕", key="add_row", help="Add new row"):
    # 添加新行时，用默认值填充
    st.session_state.rows.append(DEFAULT_ROW.copy())
    st.rerun()

# 2. 下载图标按钮（带hover提示，点击上传，还原你要的样式）
download_svg = """
<div class="download-btn">
    <svg class="download-icon" viewBox="0 0 24 24" fill="none" xmlns="http://www.w3.org/2000/svg">
        <path d="M21 15V19C21 19.5304 20.7893 20.0391 20.4142 20.4142C20.0391 20.7893 19.5304 21 19 21H5C4.46957 21 3.96086 20.7893 3.58579 20.4142C3.21071 20.0391 3 19.5304 3 19V15M17 8L12 13L7 8M12 13V3" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round"/>
    </svg>
</div>
"""
st.markdown(download_svg, unsafe_allow_html=True)

# 点击图标触发上传
if st.button("", key="upload_btn", label_visibility="collapsed"):
    st.session_state.show_upload = True
    st.rerun()

st.markdown("</div></div>", unsafe_allow_html=True)

# ====================== 6. 隐藏式上传 ======================
if st.session_state.show_upload:
    uploaded_file = st.file_uploader("Upload CSV", type="csv", label_visibility="collapsed", key="csv_upload")
    if uploaded_file is not None:
        df = pd.read_csv(uploaded_file)
        df.columns = ["spacer", "scaffold", "template", "pbs", "linker", "motif"]
        # 自动填入对应行
        for i, row in df.iterrows():
            if i < len(st.session_state.rows):
                st.session_state.rows[i]["spacer"] = str(row["spacer"])
                st.session_state.rows[i]["scaffold"] = str(row["scaffold"])
                st.session_state.rows[i]["template"] = str(row["template"])
                st.session_state.rows[i]["pbs"] = str(row["pbs"])
                st.session_state.rows[i]["linker"] = str(row["linker"])
                st.session_state.rows[i]["motif"] = str(row["motif"])
        st.session_state.show_upload = False
        st.rerun()

# ====================== 7. START按钮+运行中提示（修复iloc错误） ======================
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
                # 提取最优结果，写入Linker（list取索引，无需iloc）
                updated_rows[i]["linker"] = result[0]['linker']
            
            # 更新会话状态，刷新界面
            st.session_state.rows = updated_rows
            st.success("✅ Calculation completed!")
            st.rerun()
        except Exception as e:
            st.error(f"❌ Error: {str(e)}")

st.markdown("</div>", unsafe_allow_html=True)

# ====================== 8. 同步会话状态 ======================
st.session_state.rows = updated_rows
