import streamlit as st
st.set_page_config(page_title="pegLIT", layout="wide")

import pandas as pd
import RNA
import peglit_min

# ====================== 1. 只在第一次加载时初始化 ======================
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

# ====================== 样式 ======================
st.markdown("""
<style>
.action-row {
    display: flex;
    align-items: center;
    gap: 14px;
    margin: 8px 0;
}
.upload-icon-btn {
    width: 44px;
    height: 44px;
    display: flex;
    align-items: center;
    justify-content: center;
    border-radius: 10px;
    cursor: pointer;
}
.upload-icon-btn:hover {
    background: #f1f3f4;
}
div[data-testid="stFileUploader"] {
    display: none !important;
}
#MainMenu, footer, header {
    visibility: hidden;
}
</style>
""", unsafe_allow_html=True)

# ====================== 标题 ======================
st.markdown("# pegLIT")
st.caption("Automatically identify non-interfering nucleotide linkers between a pegRNA and 3' motif.")
st.divider()

# ====================== 表头 ======================
col_w = [1, 2, 3, 1, 1, 2]
c = st.columns(col_w)
c[0].markdown("**Spacer**")
c[1].markdown("**Scaffold**")
c[2].markdown("**Template**")
c[3].markdown("**PBS**")
c[4].markdown("**Linker**")
c[5].markdown("**Motif**")

# ====================== 表格：只读取，不覆盖 ======================
for i, row in enumerate(st.session_state.rows):
    c = st.columns(col_w)
    c[0].text_input("Spacer",    row["spacer"],    key=f"sp_{i}")
    c[1].text_input("Scaffold",  row["scaffold"],  key=f"sc_{i}")
    c[2].text_input("Template",  row["template"],  key=f"te_{i}")
    c[3].text_input("PBS",       row["pbs"],       key=f"pb_{i}")
    c[4].text_input("Linker",    row["linker"],    disabled=True, key=f"li_{i}")
    c[5].text_input("Motif",     row["motif"],     key=f"mo_{i}")

# ====================== 按钮：加号 + 上传 ======================
st.markdown("<div class='action-row'>", unsafe_allow_html=True)

# 加号
if st.button("⊕", key="add_row"):
    st.session_state.rows.append({
        "spacer": "",
        "scaffold": "",
        "template": "",
        "pbs": "",
        "linker": "NNNNNNNN",
        "motif": ""
    })
    st.rerun()

# 上传图标
st.markdown("""
<div class="upload-icon-btn" title="Import CSV">
<svg width="22" height="22" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2">
  <path d="M21 15v4a2 2 0 0 1-2 2H5a2 2 0 0 1-2-2v-4"></path>
  <path d="M7 10l5-5 5 5"></path>
  <path d="M12 15V3"></path>
</svg>
</div>
""", unsafe_allow_html=True)

# 点击触发上传
if st.button(" ", key="trigger_upload", label_visibility="collapsed"):
    st.session_state.show_upload = True
    st.rerun()

st.markdown("</div>", unsafe_allow_html=True)

# ====================== 上传CSV ======================
if st.session_state.show_upload:
    f = st.file_uploader("Upload CSV", type="csv", label_visibility="collapsed")
    if f:
        df = pd.read_csv(f)
        for i, r in df.iterrows():
            if i < len(st.session_state.rows):
                st.session_state.rows[i]["spacer"]   = str(r.iloc[0])
                st.session_state.rows[i]["scaffold"] = str(r.iloc[1])
                st.session_state.rows[i]["template"] = str(r.iloc[2])
                st.session_state.rows[i]["pbs"]      = str(r.iloc[3])
                st.session_state.rows[i]["linker"]   = str(r.iloc[4])
                st.session_state.rows[i]["motif"]    = str(r.iloc[5])
        st.session_state.show_upload = False
        st.rerun()

st.divider()

# ====================== START 运行 ======================
if st.button("START", type="primary", use_container_width=True):
    with st.spinner("Running..."):
        try:
            for i, row in enumerate(st.session_state.rows):
                res = peglit_min.pegLIT(
                    seq_spacer=row["spacer"],
                    seq_scaffold=row["scaffold"],
                    seq_template=row["template"],
                    seq_pbs=row["pbs"],
                    seq_motif=row["motif"],
                    linker_pattern=row["linker"]
                )
                # 把结果写回 session_state
                st.session_state.rows[i]["linker"] = res[0]["linker"]
            st.success("Calculation completed!")
            st.rerun()
        except Exception as e:
            st.error(f"Error: {str(e)}")
