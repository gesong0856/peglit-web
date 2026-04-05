import streamlit as st
st.set_page_config(page_title="pegLIT", layout="wide")
# 强制清除所有缓存，彻底杜绝缓存问题
st.cache_data.clear()
st.cache_resource.clear()

import pandas as pd
import RNA
import peglit_min
import random
import time

# ====================== 1. 初始化（增强健壮性 + 状态锁） ======================
DEFAULT_SEQ = {
    "spacer": "",
    "scaffold": "GTTTCAGAGCTATGCTGGAAACAGCATAGCAAGTTGAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGC",
    "template": "",
    "pbs": "",
    "linker": "NNNNNNNN",
    "motif": ""
}

if "rows" not in st.session_state:
    st.session_state.rows = [DEFAULT_SEQ.copy()]

# 🎯 核心状态锁：控制只有手动点击START才计算
if "run_calculation" not in st.session_state:
    st.session_state.run_calculation = False

# ====================== 2. 全局样式 ======================
st.markdown("""
<style>
* {
    font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif !important;
    box-sizing: border-box;
    margin: 0;
    padding: 0;
}
body {
    background-color: #ffffff;
}
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
    color: #1e40af;
    margin-bottom: 2rem;
    line-height: 1.6;
}
.table-card {
    max-width: 1200px;
    margin: 0 auto 1rem;
    border: 1px solid #e5e7eb;
    border-radius: 8px;
    overflow: hidden;
    background: white;
    box-shadow: 0 2px 8px rgba(0,0,0,0.1);
}
.table-header {
    display: grid;
    grid-template-columns: 1fr 1.5fr 1.5fr 1fr 1fr 1.5fr;
    gap: 0;
    background: #fff;
    padding: 1.25rem 0;
    font-weight: 500;
    font-size: 1.2rem;
    border-bottom: 1px solid #e5e7eb;
}
.table-header > div {
    padding: 0 1rem;
    text-align: left;
    line-height: 1.5;
}
.table-input-row {
    display: grid;
    grid-template-columns: 1fr 1.5fr 1.5fr 1fr 1fr 1.5fr;
    gap: 0;
    border-bottom: 1px solid #e5e7eb;
    align-items: center;
}
.table-input-row input {
    width: 100%;
    border: none !important;
    outline: none !important;
    font-size: 1rem;
    padding: 0.75rem 1rem;
    background: transparent !important;
    line-height: 1.5;
}
.table-input-row input:disabled {
    background: #f3f4f6 !important;
    color: #1e40af !important;
    cursor: not-allowed !important;
    opacity: 1 !important;
}
.action-row {
    display: flex;
    align-items: center;
    gap: 12px;
    padding: 12px 16px;
    height: 48px;
}
.stButton>button[kind="primary"] {
    background-color: #2563eb;
    color: white !important;
    border: none !important;
    border-radius: 12px !important;
    padding: 0.8rem 2.5rem !important;
    font-size: 1.2rem !important;
}
.stButton>button[kind="primary"]:hover {
    background-color: #1d4ed8;
}
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
<div class="table-header">
    <div>Spacer</div>
    <div>Scaffold</div>
    <div>Template</div>
    <div>PBS</div>
    <div>Linker Pattern</div>
    <div>Motif</div>
</div>
""", unsafe_allow_html=True)

for idx, row in enumerate(st.session_state.rows):
    st.markdown("<div class='table-input-row'>", unsafe_allow_html=True)
    cols = st.columns([1, 1.5, 1.5, 1, 1, 1.5])
    cols[0].text_input(label=f"spacer_{idx}", value=row["spacer"], label_visibility="collapsed", key=f"spacer_{idx}")
    cols[1].text_input(label=f"scaffold_{idx}", value=row["scaffold"], label_visibility="collapsed", key=f"scaffold_{idx}")
    cols[2].text_input(label=f"template_{idx}", value=row["template"], label_visibility="collapsed", key=f"template_{idx}")
    cols[3].text_input(label=f"pbs_{idx}", value=row["pbs"], label_visibility="collapsed", key=f"pbs_{idx}")
    cols[4].text_input(label=f"linker_{idx}", value=row["linker"], label_visibility="collapsed", disabled=True, key=f"linker_{idx}")
    cols[5].text_input(label=f"motif_{idx}", value=row["motif"], label_visibility="collapsed", key=f"motif_{idx}")
    st.markdown("</div>", unsafe_allow_html=True)

# ====================== 5. 操作按钮 ======================
st.markdown("<div class='action-row'>", unsafe_allow_html=True)
if st.button("⊕", key="add_row"):
    st.session_state.rows.append(DEFAULT_SEQ.copy())
    st.rerun()

uploaded_file = st.file_uploader("Upload CSV", type="csv", label_visibility="collapsed", key="csv_upload")
if uploaded_file is not None:
    try:
        df = pd.read_csv(uploaded_file)
        df.columns = df.columns.str.lower()
        required_cols = ["spacer", "scaffold", "template", "pbs", "linker", "motif"]
        for c in required_cols:
            if c not in df.columns:
                df[c] = DEFAULT_SEQ[c]
        df = df[required_cols]
        for i, (_, r) in enumerate(df.iterrows()):
            d = r.to_dict()
            if i < len(st.session_state.rows):
                st.session_state.rows[i] = d
            else:
                st.session_state.rows.append(d)
        st.success("✅ CSV imported")
        st.rerun()
    except Exception as e:
        st.error(f"Import failed: {e}")
st.markdown("</div></div>", unsafe_allow_html=True)

# ====================== 6. ✅ 100% 修复版 START 按钮 ======================
st.markdown("<div class='start-btn-container'>", unsafe_allow_html=True)

# 点击按钮：仅设置状态，不计算
if st.button("START", type="primary"):
    # 1. 先清空所有旧结果，彻底杜绝秒显
    for i in range(len(st.session_state.rows)):
        st.session_state.rows[i]["linker"] = "NNNNNNNN"
    # 2. 设置计算状态为True
    st.session_state.run_calculation = True
    # 3. 强制刷新页面，确保旧结果消失
    st.rerun()

# 只有状态为True时，才执行真正的计算
if st.session_state.run_calculation:
    with st.spinner("🔄 Running pegLIT... Please wait (this may take 10-30 seconds)"):
        try:
            # 强制等待，确保页面加载完成，避免缓存
            time.sleep(1)
            
            for i, r in enumerate(st.session_state.rows):
                # 🔧 严格序列预处理：只保留ATCG，大写，去空格，和官网完全一致
                def clean_seq(s):
                    return "".join([c.upper() for c in s.strip() if c.upper() in "ATCG"])

                spacer = clean_seq(r.get("spacer", ""))
                scaffold = clean_seq(r.get("scaffold", ""))
                template = clean_seq(r.get("template", ""))
                pbs = clean_seq(r.get("pbs", ""))
                motif = clean_seq(r.get("motif", ""))
                linker_pattern = "NNNNNNNN"

                st.write(f"Calculating Row {i+1}...")

                # 🎯 核心：100% 对齐官网参数，不做任何修改
                # 官网源码默认参数完全一致，不手动传seed，用默认的2020
                result = peglit_min.pegLIT(
                    seq_spacer=spacer,
                    seq_scaffold=scaffold,
                    seq_template=template,
                    seq_pbs=pbs,
                    seq_motif=motif,
                    linker_pattern=linker_pattern,
                    ac_thresh=0.5,
                    u_thresh=3,
                    n_thresh=3,
                    topn=100,
                    epsilon=1e-2,
                    num_repeats=10,
                    num_steps=250,
                    temp_init=0.15,
                    temp_decay=0.95,
                    bottleneck=1
                )

                st.write(f"Raw algorithm result: {result}")

                # 🔧 结果解析：完全复刻官网逻辑
                new_linker = "NNNNNNNN"
                if isinstance(result, list) and len(result) > 0:
                    # 官网取top1最优结果
                    if isinstance(result[0], dict):
                        new_linker = result[0].get("linker", "NNNNNNNN")
                    else:
                        new_linker = result[0]
                elif isinstance(result, str):
                    new_linker = result
                elif isinstance(result, dict):
                    new_linker = result.get("linker", "NNNNNNNN")

                # 更新结果
                st.session_state.rows[i]["linker"] = new_linker

                if new_linker == "NNNNNNNN":
                    st.warning(f"Row {i+1}: No valid linker result")
                else:
                    st.success(f"Row {i+1}: Final linker → {new_linker}")

            # 计算完成，重置状态锁
            st.session_state.run_calculation = False
            st.success("✅ Calculation completed!")
            st.rerun()

        except Exception as e:
            st.error(f"❌ Calculation error: {str(e)}")
            # 出错也重置状态
            st.session_state.run_calculation = False
            st.rerun()

st.markdown("</div>", unsafe_allow_html=True)
