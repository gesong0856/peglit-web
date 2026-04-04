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

# ================== 样式：1:1复刻官网有线表格 ==================
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
.table-card {
    max-width: 1200px;
    margin: 0 auto 1rem;
    border: 1px solid #e5e7eb;
    border-radius: 10px;
    overflow: hidden;
    background: white;
    box-shadow: 0 2px 8px rgba(0,0,0,0.1);
}

/* 表头行：有线分割+完全对齐 */
.table-header {
    display: grid;
    grid-template-columns: 1fr 1.5fr 1.5fr 1fr 1fr 1.5fr;
    background-color: #ffffff;
    padding: 1.25rem 1rem;
    font-weight: 500;
    font-size: 1.25rem;
    border-bottom: 1px solid #e5e7eb;
    /* 给表头加竖线分割 */
    background-image: linear-gradient(to right, #e5e7eb 1px, transparent 1px);
    background-size: calc(100% / 6) 100%;
    background-repeat: repeat-x;
}

/* 输入行：有线分割+完全对齐 */
.table-input-row {
    display: grid;
    grid-template-columns: 1fr 1.5fr 1.5fr 1fr 1fr 1.5fr;
    padding: 0.75rem 1rem;
    border-bottom: 1px solid #e5e7eb;
    align-items: center;
    /* 给输入行加竖线分割，和表头完全一致 */
    background-image: linear-gradient(to right, #e5e7eb 1px, transparent 1px);
    background-size: calc(100% / 6) 100%;
    background-repeat: repeat-x;
}

/* 输入框样式：官网透明无边框，适配有线表格 */
.table-input-row input {
    width: 100%;
    border: none;
    outline: none;
    font-size: 1.1rem;
    padding: 0.5rem;
    background-color: transparent;
}
.table-input-row input:focus {
    background-color: #f8f9fa;
    border-radius: 4px;
}

/* 操作按钮行：圆圈加号 + 上传箭头 */
.action-row {
    display: grid;
    grid-template-columns: 1fr 1.5fr 1.5fr 1fr 1fr 1.5fr;
    padding: 0.75rem 1rem;
    align-items: center;
}

/* 圆圈加号按钮 */
.circle-btn {
    width: 36px;
    height: 36px;
    border-radius: 50%;
    border: 2px solid #d1d5db;
    background: white;
    font-size: 20px;
    color: #6b7280;
    cursor: pointer;
    display: flex;
    align-items: center;
    justify-content: center;
}
.circle-btn:hover {
    border-color: #415E9B;
    color: #415E9B;
}

/* 上传区样式（官网拖拽区） */
.upload-area {
    max-width: 1200px;
    margin: 0 auto 2rem;
    border: 2px dashed #e5e7eb;
    border-radius: 10px;
    padding: 1.5rem;
    background-color: #f9fafb;
    display: flex;
    align-items: center;
    gap: 1rem;
}
.upload-icon {
    font-size: 1.5rem !important;
    color: #9ca3af;
    line-height: 1;
}
.upload-text {
    flex: 1;
}
.upload-text h3 {
    margin: 0 0 0.25rem;
    font-size: 1.25rem;
    color: #374151;
}
.upload-text p {
    margin: 0;
    font-size: 0.9rem;
    color: #6b7280;
}

/* START按钮（官网同款） */
.start-btn-container {
    text-align: center;
    margin: 1rem 0 3rem 0;
}
.stButton>button[kind="primary"] {
    background-color: #415E9B !important;
    color: white !important;
    border: none !important;
    border-radius: 8px !important;
    padding: 1rem 3rem !important;
    font-size: 1.75rem !important;
    font-weight: 600 !important;
    box-shadow: 0 2px 8px rgba(0,0,0,0.15) !important;
    margin: 0 auto !important;
    display: block !important;
    width: auto !important;
}
.stButton>button[kind="primary"]:hover {
    background-color: #324b7a !important;
}

/* 隐藏Streamlit默认元素 */
#MainMenu, footer, header {visibility: hidden;}
</style>
""", unsafe_allow_html=True)

# ================== 页面标题（官网一模一样） ==================
st.markdown("<h1>pegLIT</h1>", unsafe_allow_html=True)
st.markdown("""
<div class="subtitle">
Automatically identify non-interfering nucleotide<br>
linkers between a pegRNA and 3' motif.
<br><br>
<a href="#">Learn more...</a>
</div>
""", unsafe_allow_html=True)

# ================== 有线分割表格（官网布局） ==================
st.markdown("<div class='table-card'>", unsafe_allow_html=True)

# 第一行：表头
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

# 第二行：输入框
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

# 第三行：操作按钮（圆圈加号 + 上传箭头）
st.markdown("<div class='action-row'>", unsafe_allow_html=True)
col_add, col_upload, _, _, _, _ = st.columns([1, 1.5, 1.5, 1, 1, 1.5])
with col_add:
    if st.button("⊕", key="add_row", help="Add row"):
        st.session_state.rows.append({
            "spacer": "", "scaffold": "", "template": "", "pbs": "", "linker": "NNNNNNNN", "motif": ""
        })
        st.rerun()
with col_upload:
    uploaded_file = st.file_uploader("⬆️ Import CSV", type="csv", label_visibility="collapsed")
    if uploaded_file is not None:
        df = pd.read_csv(uploaded_file)
        df.columns = ["spacer", "scaffold", "template", "pbs", "linker", "motif"]
        st.session_state.rows = df.to_dict("records")
        st.rerun()
st.markdown("</div></div>", unsafe_allow_html=True)

# ================== START按钮（表格下方） ==================
st.markdown("<div class='start-btn-container'>", unsafe_allow_html=True)
if st.button("START", type="primary"):
    # 输入校验
    df_input = pd.DataFrame(updated_rows)
    if df_input.isnull().values.any() or (df_input == "").values.any():
        st.error("❌ Please fill in all fields!")
        st.stop()

    with st.spinner("🔬 Calculating..."):
        try:
            # 批量计算，并把结果直接写回 linker 框
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
                # 把最优结果直接写入 Linker Pattern
                best_linker = result.iloc[0]['linker']
                updated_rows[i]["linker"] = best_linker

            # 保存回界面
            st.session_state.rows = updated_rows
            st.success("✅ Calculation finished! Result is in Linker Pattern box.")
            st.rerun()

        except Exception as e:
            st.error(f"❌ Error: {str(e)}")
            st.exception(e)

st.markdown("</div>", unsafe_allow_html=True)

# ================== 上传区（官网样式） ==================
st.markdown("""
<div class="upload-area">
    <div class="upload-icon">☁️⬆️</div>
    <div class="upload-text">
        <h3>Drag and drop file here</h3>
        <p>Limit 200MB per file • CSV</p>
    </div>
</div>
""", unsafe_allow_html=True)

# 更新会话状态
st.session_state.rows = updated_rows
