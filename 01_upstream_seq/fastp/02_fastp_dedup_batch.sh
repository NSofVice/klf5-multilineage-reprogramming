#!/bin/bash

# ==============================================
# fastp 批量去PCR重复 + 质控脚本（并行版本）
# 功能：对指定目录下所有PE fastq.gz文件进行去重和质控
# 说明：所有输出文件均保存在当前目录下
# 日期：2026-05-14
# ==============================================

# --- 配置参数 ---
INPUT_DIR="/data/input/Files/WXY/7CL_CUTTag_clean"  # 只读目录
OUTPUT_DIR="./TPRX1_fastp_dedup_results"                   # 当前目录下的输出文件夹
THREADS=10                                                # 每个fastp任务使用的线程数
MAX_JOBS=6                                               # 同时运行的最大任务数

# 创建输出目录
mkdir -p ${OUTPUT_DIR}

# 日志文件
LOG_FILE="./fastp_batch.log"
echo "任务开始时间: $(date)" > ${LOG_FILE}

# --- fastp 参数说明 ---
# -D / --dedup                    启用去PCR重复功能
# -w / --thread                   线程数
# -j / --json                     输出JSON格式报告
# -h / --html                     输出HTML格式报告
# -l / --length_required          最小read长度，默认15
# -q / --qualified_quality_phred  高质量碱基阈值，默认15
# -u / --unqualified_percent_limit 允许的低质量碱基最大比例，默认40
# -n / --n_base_limit             N碱基最大数量，默认5

# --- 处理函数 ---
process_sample() {
    local fq1=$1
    local fq2=$2
    local sample_name=$3
    local job_id=$4
    
    echo "[$(date)] [Job-${job_id}] 开始处理样本: ${sample_name}" | tee -a ${LOG_FILE}
    
    # 定义输出文件路径（全部放在 OUTPUT_DIR 下）
    local out1="${OUTPUT_DIR}/${sample_name}_1.dedup.fq.gz"
    local out2="${OUTPUT_DIR}/${sample_name}_2.dedup.fq.gz"
    local json_report="${OUTPUT_DIR}/${sample_name}_report.json"
    local html_report="${OUTPUT_DIR}/${sample_name}_report.html"
    
    # 运行 fastp
    fastp \
        -i "${fq1}" \
        -I "${fq2}" \
        -o "${out1}" \
        -O "${out2}" \
        -D \
        -w ${THREADS} \
        -j "${json_report}" \
        -h "${html_report}" \
        -l 30 \
        -q 15 \
        -u 40 \
        -n 5 \
        --detect_adapter_for_pe \
        2>&1 | tee -a "${OUTPUT_DIR}/${sample_name}.fastp.log"
    
    # 检查运行状态
    if [ ${PIPESTATUS[0]} -eq 0 ]; then
        echo "[$(date)] [Job-${job_id}] 样本 ${sample_name} 处理成功." | tee -a ${LOG_FILE}
    else
        echo "[$(date)] [Job-${job_id}] [错误] 样本 ${sample_name} 处理失败！" | tee -a ${LOG_FILE}
        return 1
    fi
}

# 导出函数和环境变量供子shell使用
export -f process_sample
export INPUT_DIR OUTPUT_DIR THREADS LOG_FILE

# --- 主流程：查找并配对文件 ---
echo "正在扫描 ${INPUT_DIR} 中的FASTQ文件..." | tee -a ${LOG_FILE}

# 创建临时文件存储任务列表
TASK_FILE=$(mktemp)
trap "rm -f ${TASK_FILE}" EXIT

# 生成任务列表
job_id=0
while IFS= read -r -d '' fq1_file; do
    # 推断对应的 _2.clean.fq.gz 文件路径
    fq2_file="${fq1_file/_1.clean.fq.gz/_2.clean.fq.gz}"
    
    if [ -f "${fq2_file}" ]; then
        # 提取样本名
        sample_name=$(basename "${fq1_file}" | sed 's/_1\.clean\.fq\.gz$//')
        job_id=$((job_id + 1))
        echo "${fq1_file}|${fq2_file}|${sample_name}|${job_id}" >> ${TASK_FILE}
    else
        echo "[警告] 未找到配对文件: ${fq2_file}，跳过: ${fq1_file}" | tee -a ${LOG_FILE}
    fi
done < <(find ${INPUT_DIR} -type f -name "*_1.clean.fq.gz" -print0)

total_jobs=${job_id}
echo "共发现 ${total_jobs} 个样本待处理，将使用 ${MAX_JOBS} 个并行任务" | tee -a ${LOG_FILE}

# 并行执行所有任务
cat ${TASK_FILE} | xargs -I {} -P ${MAX_JOBS} bash -c '
    IFS="|" read -r fq1 fq2 sample_name job_id <<< "{}"
    process_sample "${fq1}" "${fq2}" "${sample_name}" "${job_id}"
'

echo "任务结束时间: $(date)" | tee -a ${LOG_FILE}
echo "所有任务完成，结果保存在: ${OUTPUT_DIR}" | tee -a ${LOG_FILE}

# 生成汇总报告
echo "=========================================" >> ${LOG_FILE}
echo "任务汇总（并行执行）:" >> ${LOG_FILE}
echo "  总任务数: ${total_jobs}" >> ${LOG_FILE}
echo "  并行数: ${MAX_JOBS}" >> ${LOG_FILE}
echo "  每任务线程: ${THREADS}" >> ${LOG_FILE}
echo "  总理论最大线程: $((MAX_JOBS * THREADS))" >> ${LOG_FILE}