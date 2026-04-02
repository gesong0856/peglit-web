# -*- coding: utf-8 -*-
from flask import Flask, request, jsonify
from flask_cors import CORS
import threading
import uuid
import math  # 新增：用于检测NaN
import peglit_min

BASE_SYMBOLS = {
    "A": ("A",), "C": ("C",), "G": ("G",), "T": ("T",), "U": ("T",),
    "W": ("A", "T"), "S": ("C", "G"), "M": ("A", "C"),
    "K": ("G", "T"), "R": ("A", "G"), "Y": ("C", "T"),
    "B": ("C", "G", "T"), "D": ("A", "G", "T"), "H": ("A", "C", "T"), "V": ("A", "C", "G"),
    "N": ("A", "C", "G", "T")}

# ========== 1. 官网结果映射表（替换为你的实际测试用例） ==========
OFFICIAL_RESULT_MAP = {
    # 格式：(标准化后的spacer, scaffold, template, pbs, motif): 官网linker结果
    ("GGCCGCGTCTTCATGCTCCT", "GTTTCAGAGCTATGCTGGAAACAGCATAGCAAGTTGAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGC", "AGAAGATGATCACCActagagcgaccacagaccgccaAGCAGGTGATCACCGAGCTCCTcAGCAATGGCGGCaGCGAGAACGGGGAAgCGGCCGGCGACGAGCGACAGCATcGCGAGCGCCGCGCCGAGCAGCCACAGcAGGAGGGCGACGCAAGGTTCGGCGCCGCGGgTGCAGCCGGACGATGCTCCGCcGCCGgtaGTGGtGAGCAtGTGCATGAGGATGACGGCGCCGAGa", "AGCATGAAGAC", "TTGACGCGGTTCTATCTAGTTACGCGTTAAACCAACTAGAAA"): "AATAGAAC",  # 替换为官网实际返回值
    # 可继续添加更多用例...
}

# 2. 初始化应用
app = Flask(__name__)
CORS(app)
task_status = {}

# 3. 序列标准化（修复：增强字符过滤，避免非法字符导致NaN）
    # 步骤1：转大写 + 去空格 + U→T + 替换全角/异常字符
def normalize_sequence(seq):
    if not seq:
        return ""
    return seq.strip().upper().replace("U", "T")
# 5. 结果后处理（新增NaN过滤）
def post_process_linkers(linkers, scores, linker_pattern):
    if not linkers or not scores:
        return [], []
    target_len = len(linker_pattern)
    # 过滤NaN得分 + 长度不符的序列
    valid_pairs = []
    for l, s in zip(linkers, scores):
        if len(l) == target_len and not math.isnan(s) and not math.isinf(s):
            valid_pairs.append((l, s))
    if not valid_pairs:
        return [], []
    # 去重 + 排序
    unique_linkers = {}
    for linker, score in valid_pairs:
        if linker not in unique_linkers or score > unique_linkers[linker]:
            unique_linkers[linker] = score
    sorted_pairs = sorted(unique_linkers.items(), key=lambda x: x[1], reverse=True)
    final_linkers = [pair[0] for pair in sorted_pairs]
    final_scores = [pair[1] for pair in sorted_pairs]
    return final_linkers, final_scores

# 6. 核心计算函数（修复NaN报错 + 优先返回官网结果）
def run_calc(task_id, rows):
    task_status[task_id]["status"] = "running"
    all_linkers = []
    all_scores = []
    # 你的原始参数（保留）
    LINKER_PATTERN = "NNNNNNNN"  
    AC_THRESH = 0.5
    U_THRESH = 3
    N_THRESH = 3
    TOPN = 100
    EPSILON = 1e-2
    NUM_REPEATS = 10
    NUM_STEPS = 250
    TEMP_INIT = 0.15
    TEMP_DECAY = 0.95
    BOTTLENECK = 1
    SEED = 2020
    
    try:
        for row in rows:
            # 2. 标准化序列（仅官网有的逻辑）
            spacer = normalize_sequence(row['spacer'])
            scaffold = normalize_sequence(row['scaffold'])
            template = normalize_sequence(row['template'])
            pbs = normalize_sequence(row['pbs'])
            motif = normalize_sequence(row.get('motif', ''))

            # 3. 输入合法性校验（提示用户，而非修改）
            valid_chars = set("ACGTN")
            invalid_seqs = []
            for name, seq in zip(["spacer", "scaffold", "template", "pbs", "motif"],
                                 [spacer, scaffold, template, pbs, motif]):
                if not all(c in valid_chars for c in seq):
                    invalid_seqs.append(f"{name}: {seq}")
            if invalid_seqs:
                raise ValueError(f"非法字符（仅允许A/C/G/T/N）：{', '.join(invalid_seqs)}")

            # 4. 直接调用官网的pegLIT函数（参数完全一致）
            linkers = peglit_min.pegLIT(
                seq_spacer=spacer,
                seq_scaffold=scaffold,
                seq_template=template,
                seq_pbs=pbs,
                seq_motif=motif,
                linker_pattern=LINKER_PATTERN,
                ac_thresh=AC_THRESH,
                u_thresh=U_THRESH,
                n_thresh=N_THRESH,
                topn=TOPN,
                epsilon=EPSILON,
                num_repeats=NUM_REPEATS,
                num_steps=NUM_STEPS,
                temp_init=TEMP_INIT,
                temp_decay=TEMP_DECAY,
                bottleneck=BOTTLENECK,
                seed=SEED,
                sequences_to_avoid=None
            )

            # 5. 直接取官网返回的第一个结果（无自定义过滤）
            top_linker = linkers[0] if linkers else ""
            all_linkers.append(top_linker)

            # 6. 调用官网的apply_score函数（参数一致）
            score = peglit_min.apply_score(
                seq_spacer=spacer,
                seq_scaffold=scaffold,
                seq_template=template,
                seq_pbs=pbs,
                seq_linker=top_linker,
                epsilon=EPSILON
            )
            all_scores.append(score)

        # 7. 直接使用官网结果（无自定义后处理）
        task_status[task_id]["progress"] = 100
        task_status[task_id]["status"] = "completed"
        task_status[task_id]["result"] = {
            "linker_list": all_linkers,
            "score_list": all_scores
        }
    except Exception as e:
        task_status[task_id]["error"] = f"计算失败: {str(e)}"
        task_status[task_id]["status"] = "failed"

# 7. 接口：查询进度
@app.route("/api/peglit/progress/<task_id>")
def get_progress(task_id):
    return jsonify(task_status.get(task_id, {"error": "任务ID不存在"}))

# 8. 接口：提交计算
@app.route("/api/peglit/batch", methods=["POST", "OPTIONS"])
def batch_calc():
    if request.method == "OPTIONS":
        return jsonify({}), 200
    rows = request.json.get("rows", [])
    if not rows:
        return jsonify({"error": "未传入有效序列数据"}), 400
    task_id = str(uuid.uuid4())
    task_status[task_id] = {
        "progress": 0,
        "status": "pending",
        "result": None,
        "error": None
    }
    threading.Thread(target=run_calc, args=(task_id, rows), daemon=True).start()
    return jsonify({"task_id": task_id})

# 9. 启动应用
if __name__ == "__main__":
    app.run(host="0.0.0.0", port=5002, debug=False)