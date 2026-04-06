import streamlit as st
import pandas as pd
import RNA
import peglit_min
import random
import numpy as np

# ====================== 1. 基础配置 & 初始化（修复session_state问题） ======================
st.set_page_config(page_title="pegLIT", layout="wide")

# 强制清除缓存（放在最顶部，确保生效）
st.cache_data.clear()
st.cache_resource.clear()

# 定义默认行数据
DEFAULT_SEQ = {
    "spacer": "",
    "scaffold": "GTTTCAGAGCTATGCTGGAAACAGCATAGCAAGTTGAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGC",
    "template": "",
    "pbs": "",
    "linker": "NNNNNNNN",
    "motif": "",
    "unique_seed": 2020
}

# 安全初始化session_state（避免未定义导致的错误）
if "rows" not in st.session_state:
    st.session_state.rows = [DEFAULT_SEQ.copy()]
if "run_calculation" not in st.session_state:
    st.session_state.run_calculation = False

# ====================== 2. 全局样式（完全保留） ======================
st.markdown("""
<style>
* {
    font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif !important;
    box-sizing: border-box;
    margin: 0;
    padding: 0;
}
body {background-color: #ffffff;}
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
.add-btn {
    width: 48px;
    height: 48px;
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
div[data-testid="stFileUploader"] {
    width: 48px;
    height: 48px;
    border-radius: 12px;
    border: 1px solid #d1d5db;
    background: white;
    overflow: hidden;
    position: relative;
    cursor: pointer;
    transition: all 0.2s;
}
div[data-testid="stFileUploader"]:hover {
    border-color: #3b82f6;
    background: #f3f4f6;
}
div[data-testid="stFileUploader"] > div {display: none !important;}
div[data-testid="stFileUploader"]::after {
    content: "⬆️";
    position: absolute;
    top: 50%;
    left: 50%;
    transform: translate(-50%, -50%);
    font-size: 24px;
    color: #374151;
    pointer-events: none;
}
div[data-testid="stFileUploader"]::before {
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
    pointer-events: none;
}
div[data-testid="stFileUploader"]:hover::before {
    opacity: 1;
    visibility: visible;
}
.stButton>button[kind="primary"] {
    background-color: #2563eb;
    color: white !important;
    border: none !important;
    border-radius: 12px !important;
    padding: 0.8rem 2.5rem !important;
    font-size: 1.2rem !important;
    display: block;
    margin: 0 auto;
}
.stButton>button[kind="primary"]:hover {
    background-color: #1d4ed8;
    color: white !important;
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

# 渲染输入行（增加key的唯一性，避免Streamlit交互冲突）
for idx, row in enumerate(st.session_state.rows):
    st.markdown(f"<div class='table-input-row' id='row_{idx}'>", unsafe_allow_html=True)
    cols = st.columns([1, 1.5, 1.5, 1, 1, 1.5])
    
    cols[0].text_input(
        label=f"spacer_{idx}",
        value=row["spacer"],
        label_visibility="collapsed",
        key=f"spacer_{idx}_input"  # 增加_input后缀，避免key冲突
    )
    cols[1].text_input(
        label=f"scaffold_{idx}",
        value=row["scaffold"],
        label_visibility="collapsed",
        key=f"scaffold_{idx}_input"
    )
    cols[2].text_input(
        label=f"template_{idx}",
        value=row["template"],
        label_visibility="collapsed",
        key=f"template_{idx}_input"
    )
    cols[3].text_input(
        label=f"pbs_{idx}",
        value=row["pbs"],
        label_visibility="collapsed",
        key=f"pbs_{idx}_input"
    )
    linker_disabled = len(str(row["linker"]).strip()) > 0 and row["linker"] != "NNNNNNNN"
    cols[4].text_input(
        label=f"linker_{idx}",
        value=row["linker"],
        label_visibility="collapsed",
        disabled=linker_disabled,
        key=f"linker_{idx}_input"
    )
    cols[5].text_input(
        label=f"motif_{idx}",
        value=row["motif"],
        label_visibility="collapsed",
        key=f"motif_{idx}_input"
    )
    st.markdown("</div>", unsafe_allow_html=True)
st.markdown("</div>", unsafe_allow_html=True)

# ====================== 5. 操作按钮（修复新增行逻辑） ======================
st.markdown("<div class='action-row'>", unsafe_allow_html=True)
if st.button("⊕ Add Row", key="add_row_btn", help="Add new row"):
    new_row = DEFAULT_SEQ.copy()
    new_row["unique_seed"] = max([r["unique_seed"] for r in st.session_state.rows]) + 100
    st.session_state.rows.append(new_row)
    st.rerun()  # 强制刷新，确保新增行立即显示

# CSV导入功能（简化逻辑，避免导入错误导致按钮无响应）
uploaded_file = st.file_uploader("Import CSV", type="csv", key="csv_uploader", label_visibility="collapsed")
if uploaded_file is not None:
    try:
        df = pd.read_csv(uploaded_file)
        df.columns = df.columns.str.lower()
        required_cols = ["spacer", "scaffold", "template", "pbs", "linker", "motif"]
        for col in required_cols:
            if col not in df.columns:
                df[col] = DEFAULT_SEQ[col]
        st.session_state.rows = []
        base_seed = 2020
        for i, (_, row) in enumerate(df.iterrows()):
            row_dict = DEFAULT_SEQ.copy()
            row_dict["spacer"] = str(row["spacer"]).upper().strip() if pd.notna(row["spacer"]) else ""
            row_dict["scaffold"] = str(row["scaffold"]).upper().strip() if pd.notna(row["scaffold"]) else DEFAULT_SEQ["scaffold"]
            row_dict["template"] = str(row["template"]).upper().strip() if pd.notna(row["template"]) else ""
            row_dict["pbs"] = str(row["pbs"]).upper().strip() if pd.notna(row["pbs"]) else ""
            row_dict["linker"] = str(row["linker"]).upper().strip() if pd.notna(row["linker"]) else "NNNNNNNN"
            row_dict["motif"] = str(row["motif"]).upper().strip() if pd.notna(row["motif"]) else ""
            row_dict["unique_seed"] = base_seed + i * 100
            st.session_state.rows.append(row_dict)
        st.success("✅ CSV imported successfully!")
        st.rerun()
    except Exception as e:
        st.error(f"❌ CSV import failed: {str(e)}")
st.markdown("</div>", unsafe_allow_html=True)

# ====================== 6. 核心工具函数（简化逻辑，减少执行错误） ======================
def validate_nucleotide_seq(seq):
    valid_chars = set("ATCGN")
    seq = seq.upper().strip()
    if not seq:
        return True, ""
    if set(seq) - valid_chars:
        return False, f"Invalid chars: {set(seq) - valid_chars} (only A/T/C/G/N allowed)"
    return True, ""

def preprocess_for_peglit(spacer, scaffold, template, pbs, motif, linker_pattern, seed=2020):
    random.seed(seed)
    np.random.seed(seed)
    # 标准化序列
    spacer = spacer.upper().strip().replace("T", "U")
    scaffold = scaffold.upper().strip().replace("T", "U")
    template = template.upper().strip().replace("T", "U")
    pbs = pbs.upper().strip().replace("T", "U")
    motif = motif.upper().strip().replace("T", "U")
    linker_pattern = linker_pattern.upper().strip().replace("T", "U") or "NNNNNNNN"
    # 计算AC阈值
    ac_thresh = 0.5 * len(linker_pattern)
    return {
        "spacer": spacer, "scaffold": scaffold, "template": template,
        "pbs": pbs, "motif": motif, "linker_pattern": linker_pattern,
        "ac_thresh": ac_thresh, "seed": seed
    }

# ====================== 7. START按钮（核心修复：独立触发逻辑+详细日志） ======================
# 单独的START按钮区域，确保触发优先级
col_start, _, _ = st.columns([1, 2, 1])
with col_start:
    if st.button("START", type="primary", key="start_calc_btn"):
        # 标记开始计算，避免Streamlit重渲染导致的重复执行
        st.session_state.run_calculation = True
        # 立即显示加载状态
        with st.spinner("🔄 Running calculations..."):
            try:
                total_rows = len(st.session_state.rows)
                used_linkers = set()
                success_count = 0

                # 遍历每一行计算（实时更新状态，避免卡顿）
                for i in range(total_rows):
                    st.write(f"### Processing Row {i+1}/{total_rows}")
                    # 读取当前行输入（从input key读取，确保获取最新值）
                    row = st.session_state.rows[i]
                    spacer = st.session_state[f"spacer_{i}_input"]
                    scaffold = st.session_state[f"scaffold_{i}_input"]
                    template = st.session_state[f"template_{i}_input"]
                    pbs = st.session_state[f"pbs_{i}_input"]
                    motif = st.session_state[f"motif_{i}_input"]
                    linker_pattern = st.session_state[f"linker_{i}_input"] or "NNNNNNNN"
                    unique_seed = row["unique_seed"]

                    # 1. 校验序列合法性
                    all_seq = spacer + scaffold + template + pbs + motif + linker_pattern
                    is_valid, err_msg = validate_nucleotide_seq(all_seq)
                    if not is_valid:
                        st.error(f"Row {i+1}: ❌ {err_msg}")
                        row["linker"] = f"NNNNNNNN_{i}"
                        used_linkers.add(row["linker"])
                        continue

                    # 2. 校验核心序列非空
                    if not all([spacer, scaffold, template, pbs, motif]):
                        st.warning(f"Row {i+1}: ⚠️ Core sequences (spacer/scaffold/template/pbs/motif) cannot be empty!")
                        row["linker"] = f"NNNNNNNN_{i}"
                        used_linkers.add(row["linker"])
                        continue

                    # 3. 预处理
                    preprocessed = preprocess_for_peglit(
                        spacer, scaffold, template, pbs, motif, linker_pattern, unique_seed
                    )

                    # 4. 调用pegLIT算法（核心：增加详细异常捕获）
                    base_linker = ""
                    try:
                        # 转回T传给算法（pegLIT预期DNA序列）
                        result = peglit_min.pegLIT(
                            seq_spacer=preprocessed["spacer"].replace("U", "T"),
                            seq_scaffold=preprocessed["scaffold"].replace("U", "T"),
                            seq_template=preprocessed["template"].replace("U", "T"),
                            seq_pbs=preprocessed["pbs"].replace("U", "T"),
                            seq_motif=preprocessed["motif"].replace("U", "T"),
                            linker_pattern=preprocessed["linker_pattern"].replace("U", "T"),
                            ac_thresh=preprocessed["ac_thresh"],
                            u_thresh=3, n_thresh=3, topn=100,
                            epsilon=1e-2, num_repeats=10, num_steps=250,
                            temp_init=0.15, temp_decay=0.95, bottleneck=1,
                            seed=preprocessed["seed"], sequences_to_avoid=None
                        )
                        # 解析结果
                        if isinstance(result, str):
                            base_linker = result.strip().replace("U", "T")
                        elif isinstance(result, list) and len(result) > 0:
                            base_linker = result[0].get("linker", "").strip().replace("U", "T") if isinstance(result[0], dict) else str(result[0]).strip().replace("U", "T")
                        elif isinstance(result, dict):
                            base_linker = result.get("linker", "").strip().replace("U", "T")
                    except ValueError as e:
                        if "not enough values to unpack" in str(e):
                            st.error(f"Row {i+1}: ❌ No valid linker (filter too strict)")
                            base_linker = f"NNNNNNNN_{i}"
                        else:
                            st.error(f"Row {i+1}: ❌ Value Error - {str(e)}")
                            base_linker = f"NNNNNNNN_{i}"
                    except Exception as e:
                        st.error(f"Row {i+1}: ❌ Calculation Error - {str(e)}")
                        base_linker = f"NNNNNNNN_{i}"

                    # 5. 确保linker唯一
                    final_linker = base_linker
                    retry = 0
                    while final_linker in used_linkers and retry < 3:
                        retry += 1
                        new_seed = unique_seed + retry * 10
                        random.seed(new_seed)
                        np.random.seed(new_seed)
                        # 重新计算
                        try:
                            result = peglit_min.pegLIT(
                                seq_spacer=preprocessed["spacer"].replace("U", "T"),
                                seq_scaffold=preprocessed["scaffold"].replace("U", "T"),
                                seq_template=preprocessed["template"].replace("U", "T"),
                                seq_pbs=preprocessed["pbs"].replace("U", "T"),
                                seq_motif=preprocessed["motif"].replace("U", "T"),
                                linker_pattern=preprocessed["linker_pattern"].replace("U", "T"),
                                ac_thresh=preprocessed["ac_thresh"],
                                u_thresh=3, n_thresh=3, topn=100,
                                epsilon=1e-2, num_repeats=10, num_steps=250,
                                temp_init=0.15, temp_decay=0.95, bottleneck=1,
                                seed=new_seed, sequences_to_avoid=None
                            )
                            final_linker = result.strip().replace("U", "T") if isinstance(result, str) else f"NNNNNNNN_{i}_{retry}"
                        except:
                            final_linker = f"NNNNNNNN_{i}_{retry}"

                    # 6. 更新结果
                    row["linker"] = final_linker
                    used_linkers.add(final_linker)
                    st.session_state.rows[i] = row  # 强制更新session_state

                    # 7. 结果提示
                    if not final_linker.startswith("NNNNNNNN_"):
                        st.success(f"Row {i+1}: ✅ Linker - {final_linker}")
                        success_count += 1
                    else:
                        st.warning(f"Row {i+1}: ⚠️ Default Linker - {final_linker}")

                # 最终汇总
                st.success(f"🎉 Calculations done! {success_count}/{total_rows} valid linkers generated.")
                # 重置计算标记
                st.session_state.run_calculation = False
                # 强制刷新，显示最新linker结果
                st.rerun()

            except Exception as e:
                st.error(f"❌ Global Error: {str(e)}")
                st.exception(e)  # 显示完整错误栈，方便排查
                st.session_state.run_calculation = False
