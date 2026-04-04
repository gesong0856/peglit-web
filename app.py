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

# ====================== 2. 全局样式（彻底解决错位） ======================
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
    font-size: 2.5rem;
    font-weight: 700;
    margin: 2rem 0 0.5rem !important;
    color: #1f2937;
}
.subtitle {
    text-align: center;
    font-size: 1rem;
    color: #6b7280;
    margin-bottom: 2rem;
}

/* 表格卡片 */
.table-card {
    max-width: 1400px;
    margin: 0 auto 1rem;
    border: 1px solid #e5e7eb;
    border-radius: 8px;
    background: white;
    box-shadow: 0 2px 8px rgba(0,0,0,0.05);
}

/* 表头 */
.table-header {
    display: grid;
    grid-template-columns: 1fr 1.8fr 2.5fr 1fr 1fr 1.8fr;
    padding: 1rem 0;
    font-weight: 500;
    border-bottom: 1px solid #e5e7eb;
}
.table-header > div {
    padding: 0 1rem;
    line-height: 1.5;
}

/* 输入行 */
.table-input-row {
    display: grid;
    grid-template-columns: 1fr 1.8fr 2.5fr 1fr 1fr 1.8fr;
    border-bottom: 1px solid #e5e7eb;
    align-items: center;
}
.table-input-row input {
    width: 100%;
    border: none !important;
    outline: none !important;
    padding: 0.75rem 1rem;
    background: transparent !important;
    line-height: 1.5;
}
.table-input-row input:disabled {
    background: #f3f4f6 !important;
    color: #374151 !important;
    cursor: not-allowed !important;
}

/* 按钮行（核心修复：对齐） */
.action-row {
    display: flex;
    align-items: center;
    gap: 12px;
    padding: 12px 16px;
    height: 56px; /* 统一高度，彻底解决错位 */
}

/* 加号按钮 */
.add-btn {
    width: 44px;
    height: 44px;
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

/* 上传图标按钮（与加号完美对齐） */
.upload-btn-wrapper {
    width: 44px;
    height: 44px;
    border-radius: 12px;
    display: flex;
    align-items: center;
    justify-content: center;
    cursor: pointer;
    transition: all 0.2s;
}
.upload-btn-wrapper:hover {
    background: #f3f4f4;
}
.upload-icon {
    width: 20px;
    height: 20px;
    stroke: #374151;
    stroke-width: 2;
}

/* 隐藏原生上传区 */
div[data-testid="stFileUploader"] {
    display: none !important;
}

/* START按钮 */
.start-btn-container {
    text-align: center;
    margin: 1.5rem 0;
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

# ====================== 4. 表格渲染 ======================
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
        # Linker Pattern 只读，仅输出
        "linker": cols[4].text_input(f"linker_{idx}", value=row["linker"], label_visibility="collapsed", disabled=True),
        "motif": cols[5].text_input(f"motif_{idx}", value=row["motif"], label_visibility="collapsed"),
    }
    updated_rows.append(updated_row)
    st.markdown("</div>", unsafe_allow_html=True)

# ====================== 5. 操作按钮行（加号 + 上传图标） ======================
st.markdown("<div class='action-row'>", unsafe_allow_html=True)

# 1. 加号按钮：添加新行
if st.button("⊕", key="add_row_btn", help="Add new row"):
    st.session_state.rows.append({
        "spacer": "", "scaffold": "", "template": "", "pbs": "", "linker": "NNNNNNNN", "motif": ""
    })
    st.rerun()

# 2. 上传图标按钮：点击触发上传（对齐 + 有 hover 提示）
st.markdown("""
<div class="upload-btn-wrapper" title="Import CSV">
    <svg class="upload-icon" viewBox="0 0 24 24" fill="none" xmlns="http://www.w3.org/2000/svg">
        <path d="M21 15v4a2 2 0 0 1-2 2H5a2 2 0 0 1-2-2v-4"></path>
        <path d="M7 10l5 5 5-5"></path>
        <path d="M12 15V3"></path>
    </svg>
</div>
""", unsafe_allow_html=True)

# 点击图标触发上传
if st.button("", key="upload_trigger", label_visibility="collapsed"):
    st.session_state.show_upload = True
    st.rerun()

st.markdown("</div></div>", unsafe_allow_html=True)

# ====================== 6. 隐藏式上传功能 ======================
if st.session_state.show_upload:
    uploaded_file = st.file_uploader("Upload CSV", type="csv", label_visibility="collapsed", key="csv_file")
    if uploaded_file is not None:
        df = pd.read_csv(uploaded_file)
        # 强制列名匹配
        df.columns = ["spacer", "scaffold", "template", "pbs", "linker", "motif"]
        # 自动填入对应行
        for i, row in df.iterrows():
            if i < len(st.session_state.rows):
                st.session_state.rows[i] = dict(row)
        st.session_state.show_upload = False
        st.rerun()

# ====================== 7. START 按钮 + 修复 iloc 错误 ======================
st.markdown("<div class='start-btn-container'>", unsafe_allow_html=True)
if st.button("START", type="primary"):
    with st.spinner("🔄 Running... Please wait"):
        try:
            # 遍历每一行进行计算
            for i, r in enumerate(updated_rows):
                # 调用核心函数
                result = peglit_min.pegLIT(
                    seq_spacer=r["spacer"],
                    seq_scaffold=r["scaffold"],
                    seq_template=r["template"],
                    seq_pbs=r["pbs"],
                    seq_motif=r["motif"],
                    linker_pattern=r["linker"]
                )
                # 🔥 关键修复：result 是 List，直接取索引 [0] 即可，无需 iloc
                updated_rows[i]["linker"] = result[0]['linker']
            
            # 更新会话状态并刷新
            st.session_state.rows = updated_rows
            st.success("✅ Calculation completed successfully!")
            st.rerun()
        except Exception as e:
            st.error(f"❌ Execution Error: {str(e)}")

st.session_state.rows = updated_rows
