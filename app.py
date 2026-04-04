# 【必须放在最最开头！】第一个 Streamlit 命令
import streamlit as st
st.set_page_config(page_title="pegLIT 长序列版", layout="wide")

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
        {"spacer": "", "scaffold": "", "template": "", "pbs": "", "motif": "", "linker": "NNNNNNNN"}
    ]

# ================== 样式：圆圈加号 + 箭头按钮 + START按钮 ==================
st.markdown("""
<style>
/* 圆圈加号按钮 */
.circle-btn {
    font-size: 18px !important;
    border-radius: 50% !important;
    width: 32px !important;
    height: 32px !important;
    display: inline-flex;
    align-items: center;
    justify-content: center;
    border: 1px solid #ccc;
    background: white;
    cursor: pointer;
    margin-right: 10px;
}
/* 上传图标 */
.upload-icon {
    font-size: 18px;
    margin-left: 10px;
}
/* 官网同款START按钮 */
.stButton>button[kind="primary"] {
    background-color: #415E9B !important;
    color: white !important;
    border: none !important;
    border-radius: 8px !important;
    padding: 1rem 3rem !important;
    font-size: 1.75rem !important;
    font-weight: 600 !important;
    box-shadow: 0 2px 8px rgba(0,0,0,0.15) !important;
    margin: 2rem auto !important;
    display: block !important;
    width: auto !important;
}
.stButton>button[kind="primary"]:hover {
    background-color: #324b7a !important;
}
</style>
""", unsafe_allow_html=True)

# ================== 网页界面 ==================
st.title("🧬 pegLIT 长序列工具")
st.markdown("---")

# ================== 功能按钮：圆圈加号 + 上传 CSV ==================
col_btn1, col_btn2 = st.columns([0.1, 1])
with col_btn1:
    # 圆圈加号 → 添加一行
    if st.button("⊕", help="Add row"):
        st.session_state.rows.append({
            "spacer": "", "scaffold": "", "template": "", "pbs": "", "motif": "", "linker": "NNNNNNNN"
        })
        st.rerun()

with col_btn2:
    # 箭头 → 导入 CSV
    uploaded_file = st.file_uploader("⬆️ Import CSV", type="csv", label_visibility="collapsed")
    if uploaded_file is not None:
        df = pd.read_csv(uploaded_file)
        # 确保CSV列名匹配
        df.columns = ["spacer", "scaffold", "template", "pbs", "motif", "linker"]
        st.session_state.rows = df.to_dict("records")
        st.rerun()

st.markdown("---")

# ================== 动态多行输入区域（调整输入框大小） ==================
updated_rows = []
for idx, row in enumerate(st.session_state.rows):
    st.markdown(f"#### 序列组 {idx+1}")
    col1, col2 = st.columns(2)
    
    # 把 height 从 100 调整为 60，输入框更小
    with col1:
        spacer = st.text_area(
            f"Spacer 序列", 
            value=row["spacer"], 
            height=60,  # 调小高度
            placeholder="输入Spacer序列",
            key=f"sp_{idx}"
        )
        scaffold = st.text_area(
            f"Scaffold 序列", 
            value=row["scaffold"], 
            height=60,  # 调小高度
            placeholder="输入Scaffold序列",
            key=f"sc_{idx}"
        )
        template = st.text_area(
            f"Template 序列", 
            value=row["template"], 
            height=60,  # 调小高度
            placeholder="输入Template序列",
            key=f"t_{idx}"
        )
    
    with col2:
        pbs = st.text_area(
            f"PBS 序列", 
            value=row["pbs"], 
            height=60,  # 调小高度
            placeholder="输入PBS序列",
            key=f"p_{idx}"
        )
        motif = st.text_area(
            f"Motif 序列", 
            value=row["motif"], 
            height=60,  # 调小高度
            placeholder="输入Motif序列",
            key=f"m_{idx}"
        )
        linker = st.text_input(
            f"Linker 模式", 
            value=row["linker"], 
            help="默认8个N，可自定义长度",
            key=f"lk_{idx}"
        )
    
    # 保存更新后的数据
    updated_rows.append({
        "spacer": spacer,
        "scaffold": scaffold,
        "template": template,
        "pbs": pbs,
        "motif": motif,
        "linker": linker
    })
    st.markdown("---")

# 更新会话状态（关键：避免数据丢失）
st.session_state.rows = updated_rows

# ================== 官网同款START按钮 ==================
if st.button("START", type="primary"):
    all_results = []
    # 输入校验
    for i, r in enumerate(st.session_state.rows):
        if not all([r["spacer"], r["scaffold"], r["template"], r["pbs"], r["motif"]]):
            st.error(f"❌ 第 {i+1} 组序列未填写完整！")
            st.stop()

    with st.spinner("🔬 正在批量计算，请稍候..."):
        try:
            # 批量调用pegLIT核心算法
            for r in st.session_state.rows:
                result = peglit_min.pegLIT(
                    seq_spacer=r["spacer"],
                    seq_scaffold=r["scaffold"],
                    seq_template=r["template"],
                    seq_pbs=r["pbs"],
                    seq_motif=r["motif"],
                    linker_pattern=r["linker"],
                    ac_thresh=0.5,
                    u_thresh=3,
                    n_thresh=3,
                    topn=100,
                    epsilon=1e-2,
                    num_repeats=10,
                    num_steps=250,
                    temp_init=0.15,
                    temp_decay=0.9,
                    bottleneck=1,
                    seed=2026,
                    sequences_to_avoid=None
                )
                all_results.append(result)

            # 合并结果并展示
            final_result = pd.concat(all_results, ignore_index=True)
            st.success("✅ 全部计算完成！")
            st.dataframe(final_result, use_container_width=True)

            # 下载结果
            csv = final_result.to_csv(index=False)
            st.download_button(
                label="📥 下载结果CSV",
                data=csv,
                file_name="peglit_result.csv",
                mime="text/csv",
                use_container_width=True
            )

        except Exception as e:
            st.error(f"❌ 计算出错：{str(e)}")
            st.exception(e)
