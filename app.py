# -*- coding: utf-8 -*-
# 已删除：所有 Flask / CORS / threading 无关代码
# 已删除：官网映射表、任务队列、后端接口
# 只保留：纯 Streamlit 网页 + peglit 调用

import streamlit as st
import numpy as np
import pandas as pd
import RNA
import peglit_min


# ===================== 网页界面 =====================
st.set_page_config(page_title="pegLIT 长序列版", layout="wide")
st.title("pegLIT — 支持 >150nt 长序列")
st.markdown("---")

# 输入区域
col1, col2 = st.columns(2)

with col1:
    spacer = st.text_area("Spacer 序列")
    scaffold = st.text_area("Scaffold 序列")
    template = st.text_area("Template 序列")

with col2:
    pbs = st.text_area("PBS 序列")
    motif = st.text_area("Motif 序列")
    linker_pattern = st.text_input("Linker 模式", value="NNNNNNNN")

# 按钮
if st.button("运行计算"):
    if not all([spacer, scaffold, template, pbs, motif]):
        st.error("请填写所有序列！")
        st.stop()

    with st.spinner("正在计算..."):
        try:
            # 直接调用原版 peglit
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
                bottleneck=1,
                seed=2020,
                sequences_to_avoid=None
            )

            st.success("计算完成！")
            st.dataframe(result)

            # 下载
            csv = result.to_csv(index=False)
            st.download_button("下载结果", csv, "peglit_result.csv")

        except Exception as e:
            st.error(f"错误：{str(e)}")
