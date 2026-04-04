# 【必须放在最最开头！】第一个 Streamlit 命令，不能有任何前置 st.xxx 或 import
import streamlit as st
st.set_page_config(page_title="pegLIT 长序列版", layout="wide")

# 后续导入
import numpy as np
import pandas as pd
import RNA
import peglit_min

# 新版ViennaRNA默认支持长序列，删除旧版API调用，彻底解决AttributeError
# RNA.set_parameter("max_length", 2000)

# ================== 网页界面 ==================
st.title("🧬 pegLIT 长序列工具")
st.markdown("---")

# 输入区域
col1, col2 = st.columns(2)
with col1:
    spacer = st.text_area("Spacer 序列", height=100, placeholder="输入Spacer序列")
    scaffold = st.text_area("Scaffold 序列", height=100, placeholder="输入Scaffold序列")
    template = st.text_area("Template 序列", height=100, placeholder="输入Template序列")

with col2:
    pbs = st.text_area("PBS 序列", height=100, placeholder="输入PBS序列")
    motif = st.text_area("Motif 序列", height=100, placeholder="输入Motif序列")
    linker_pattern = st.text_input("Linker 模式", value="NNNNNNNN", help="默认8个N，可自定义长度")

# 运行按钮
if st.button("🚀 运行计算", type="primary", use_container_width=True):
    # 输入校验
    if not all([spacer, scaffold, template, pbs, motif]):
        st.error("❌ 请填写所有序列！")
        st.stop()

    with st.spinner("🔬 正在计算，请稍候..."):
        try:
            # 调用pegLIT核心算法
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
                temp_decay=0.9,
                bottleneck=1,
                seed=2026,
                sequences_to_avoid=None
            )

            # 展示结果
            st.success("✅ 计算完成！")
            st.dataframe(result, use_container_width=True)

            # 下载结果
            csv = result.to_csv(index=False)
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
