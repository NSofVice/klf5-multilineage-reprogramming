#!/bin/bash
# ==============================================
# CUT&Tag TE专用分析脚本
# 特点：唯一比对(-k 1)，但保留低质量reads
# ==============================================

set -euo pipefail

# --- 配置参数 ---
INPUT_DIR="./TPRX1_fastp_dedup_results"
OUTPUT_DIR="./TRPX1_CUTTag_TE_only"
THREADS=60

MM10_INDEX="/data/input/Files/ReferenceData/bt_mm10_index/mm10_index"
ECOLI_INDEX="/data/input/Files/ReferenceData/bt_Ecoli_index/Ecoli_index"
MM10_EFFECTIVE_SIZE=2652783500

mkdir -p ${OUTPUT_DIR}/{bam,bigwig,spikein,logs}

LOG_FILE="${OUTPUT_DIR}/analysis_TE.log"
PROGRESS_FILE="${OUTPUT_DIR}/completed_samples.txt"
FAILED_SAMPLES_FILE="${OUTPUT_DIR}/failed_samples.txt"

echo "=========================================" | tee -a ${LOG_FILE}
echo "CUT&Tag TE专用分析脚本" | tee -a ${LOG_FILE}
echo "特点: -k 1唯一比对，但不过滤低质量reads" | tee -a ${LOG_FILE}
echo "任务开始时间: $(date)" | tee -a ${LOG_FILE}
echo "线程数: ${THREADS}" | tee -a ${LOG_FILE}
echo "=========================================" | tee -a ${LOG_FILE}

# ==============================================
# 函数：带超时检测的 bowtie2 运行
# ==============================================
run_bowtie2() {
    local description=$1
    shift
    local log_file=$1
    shift
    
    echo "    [$(date)] ${description} 开始..." | tee -a ${LOG_FILE}
    
    if timeout 86400 bowtie2 "$@" 2> "${log_file}"; then
        echo "    [$(date)] ${description} 完成." | tee -a ${LOG_FILE}
        return 0
    else
        local exit_code=$?
        if [[ ${exit_code} -eq 124 ]]; then
            echo "    [错误] ${description} 超时（24小时）!" | tee -a ${LOG_FILE}
        else
            echo "    [错误] ${description} 失败，退出码: ${exit_code}" | tee -a ${LOG_FILE}
        fi
        return 1
    fi
}

# ==============================================
# 函数：检查样本是否已完成
# ==============================================
is_sample_completed() {
    local sample_name=$1
    
    local required_files=(
        "${OUTPUT_DIR}/bam/${sample_name}_sorted.bam"
        "${OUTPUT_DIR}/bigwig/${sample_name}_RPKM.bw"
        "${OUTPUT_DIR}/bigwig/${sample_name}_CPM.bw"
    )
    
    for file in "${required_files[@]}"; do
        if [[ ! -s "${file}" ]]; then
            return 1
        fi
    done
    return 0
}

# ==============================================
# 函数：处理单个样本（TE专用）
# ==============================================
process_sample() {
    local fq1=$1
    local fq2=$2
    local sample_name=$3

    echo "" | tee -a ${LOG_FILE}
    echo "[$(date)] ========== 开始处理TE样本: ${sample_name} ==========" | tee -a ${LOG_FILE}

    # ==========================================
    # TE模式：唯一比对 (-k 1)，但保留低质量reads
    # ==========================================
    echo "  [TE分析] 运行 Bowtie2 唯一比对 (-k 1)..." | tee -a ${LOG_FILE}
    
    # 主基因组比对（唯一比对）
    run_bowtie2 "TE-主基因组" \
        ${OUTPUT_DIR}/logs/${sample_name}_main.mapping.log \
        -p ${THREADS} --no-mixed --no-discordant --no-unal \
        --end-to-end --very-sensitive -I 10 -X 700 --phred33 -k 1 \
        -x ${MM10_INDEX} -1 ${fq1} -2 ${fq2} \
        -S ${OUTPUT_DIR}/bam/${sample_name}_main.sam || return 1

    # Spike-in比对（同样使用-k 1）
    run_bowtie2 "TE-Spike-in" \
        ${OUTPUT_DIR}/logs/${sample_name}_spikein.mapping.log \
        -p ${THREADS} --no-mixed --no-discordant --no-unal \
        --end-to-end --very-sensitive -I 10 -X 700 --phred33 -k 1 \
        -x ${ECOLI_INDEX} -1 ${fq1} -2 ${fq2} \
        -S ${OUTPUT_DIR}/spikein/${sample_name}_spikein.sam || return 1

    # ==========================================
    # BAM处理（TE专用：不过滤低质量reads）
    # 关键：去掉 -q 30 参数
    # ==========================================
    echo "  [TE分析] BAM处理（保留所有mapped reads，不过滤质量）..." | tee -a ${LOG_FILE}
    process_bam_te "${OUTPUT_DIR}/bam/${sample_name}_main.sam" \
                   "${OUTPUT_DIR}/bam/${sample_name}" || return 1

    # ==========================================
    # 生成Spike-in统计报告
    # ==========================================
    echo "  [TE分析] 生成Spike-in统计..." | tee -a ${LOG_FILE}
    
    local spikein_rate=$(grep "overall alignment rate" "${OUTPUT_DIR}/logs/${sample_name}_spikein.mapping.log" | awk '{print $1}' | sed 's/%//')
    
    {
        echo "=== ${sample_name} TE-mode Spike-in QC Report ==="
        echo "Spike-in alignment rate: ${spikein_rate}%"
        echo ""
        echo "TE分析参数说明："
        echo "  - Bowtie2参数: -k 1 (唯一比对)"
        echo "  - 不过滤MAPQ质量分数（保留所有mapped reads）"
        echo "  - 不过滤重复reads"
        echo ""
        echo "与常规模式的区别："
        echo "  常规模式: -k 1 + 过滤MAPQ≥30"
        echo "  TE模式:   -k 1 + 不过滤低质量reads"
    } > "${OUTPUT_DIR}/spikein/${sample_name}_spikein_qc.txt"

    # ==========================================
    # 生成BigWig（TE专用：不加minMappingQuality过滤）
    # 关键：去掉 --minMappingQuality 30 参数
    # ==========================================
    echo "  [TE分析] 生成BigWig（不过滤低质量reads）..." | tee -a ${LOG_FILE}
    generate_bigwig_te "${OUTPUT_DIR}/bam/${sample_name}_sorted.bam" \
                       "${OUTPUT_DIR}/bigwig/${sample_name}" \
                       ${HG38_EFFECTIVE_SIZE} || return 1

    # ==========================================
    # 统计信息
    # ==========================================
    local total_reads=$(samtools view -c "${OUTPUT_DIR}/bam/${sample_name}_sorted.bam")
    local high_qual_reads=$(samtools view -c -q 30 "${OUTPUT_DIR}/bam/${sample_name}_sorted.bam" 2>/dev/null || echo "0")
    local low_qual_reads=$((total_reads - high_qual_reads))
    
    echo "  [TE分析] 统计信息:" | tee -a ${LOG_FILE}
    echo "    总mapped reads: ${total_reads}" | tee -a ${LOG_FILE}
    echo "    高质量reads (MAPQ≥30): ${high_qual_reads}" | tee -a ${LOG_FILE}
    echo "    低质量reads (MAPQ<30): ${low_qual_reads}" | tee -a ${LOG_FILE}
    
    echo "${sample_name}" >> "${PROGRESS_FILE}"
    echo "[$(date)] ========== TE样本 ${sample_name} 完成! ==========" | tee -a ${LOG_FILE}
    return 0
}

# ==============================================
# 函数：TE专用BAM处理
# 关键修改：去掉 -q 30 过滤
# ==============================================
process_bam_te() {
    local sam_file=$1
    local output_prefix=$2

    # 步骤1：SAM转BAM
    # 关键：只排除完全未比对的reads (-F 0x04)
    # 不过滤MAPQ分数（不加 -q 30）
    samtools view -u -bS -F 0x04 "${sam_file}" > "${output_prefix}_filtered.bam" || return 1

    # 步骤2：按名称排序（为去重做准备）
    samtools sort -n -@ ${THREADS} "${output_prefix}_filtered.bam" > "${output_prefix}_namesorted.bam" || return 1
    
    # 步骤3：去重（保留重复，不上算去除）
    # TE分析保留重复，因为重复区域对TE很重要
    # 这里跳过markdup步骤，直接使用namesorted.bam
    echo "    跳过markdup，保留所有重复reads用于TE分析" | tee -a ${LOG_FILE}
    
    # 步骤4：按位置排序
    samtools sort -@ ${THREADS} "${output_prefix}_namesorted.bam" > "${output_prefix}_sorted.bam" || return 1
    
    # 步骤5：索引
    samtools index "${output_prefix}_sorted.bam" || return 1

    # 清理中间文件
    rm -f "${sam_file}" "${output_prefix}_filtered.bam" "${output_prefix}_namesorted.bam"
    
    return 0
}

# ==============================================
# 函数：TE专用生成BigWig
# 关键修改：去掉 --minMappingQuality 30 参数
# ==============================================
generate_bigwig_te() {
    local bam_file=$1
    local output_prefix=$2
    local effective_size=$3

    # RPKM标准化（不加minMappingQuality过滤；与稿件 Methods 一致）
    bamCoverage -b ${bam_file} \
        --normalizeUsing RPKM \
        --effectiveGenomeSize ${effective_size} \
        --numberOfProcessors ${THREADS} \
        --ignoreDuplicates \
        -o "${output_prefix}_RPKM.bw" 2>/dev/null || return 1

    # CPM标准化（不加minMappingQuality过滤）
    bamCoverage -b ${bam_file} \
        --normalizeUsing CPM \
        --numberOfProcessors ${THREADS} \
        --ignoreDuplicates \
        -o "${output_prefix}_CPM.bw" 2>/dev/null || return 1

    return 0
}

# ==============================================
# 主流程
# ==============================================
shopt -s nullglob
declare -A samples

for fq1 in ${INPUT_DIR}/*_1.dedup.fq.gz; do
    fq2="${fq1/_1.dedup.fq.gz/_2.dedup.fq.gz}"
    if [[ -f "${fq2}" ]]; then
        sample_name=$(basename "${fq1}" | sed 's/_1\.dedup\.fq\.gz$//')
        samples["${sample_name}"]="${fq1}|${fq2}"
    fi
done

echo "共发现 ${#samples[@]} 个样本." | tee -a ${LOG_FILE}

touch "${PROGRESS_FILE}" "${FAILED_SAMPLES_FILE}"

for sample_name in "${!samples[@]}"; do
    mkdir -p ${OUTPUT_DIR}/logs ${OUTPUT_DIR}/spikein
    
    if is_sample_completed "${sample_name}"; then
        echo "[$(date)] 跳过已完成样本: ${sample_name}" | tee -a ${LOG_FILE}
        continue
    fi

    IFS='|' read -r fq1 fq2 <<< "${samples[${sample_name}]}"
    if ! process_sample "${fq1}" "${fq2}" "${sample_name}"; then
        echo "${sample_name}" >> "${FAILED_SAMPLES_FILE}"
        echo "[$(date)] 样本 ${sample_name} 失败，继续下一个..." | tee -a ${LOG_FILE}
    fi
done

echo "=========================================" | tee -a ${LOG_FILE}
echo "任务结束时间: $(date)" | tee -a ${LOG_FILE}
echo "已完成: $(wc -l < ${PROGRESS_FILE} 2>/dev/null || echo 0)" | tee -a ${LOG_FILE}
echo "失败: $(wc -l < ${FAILED_SAMPLES_FILE} 2>/dev/null || echo 0)" | tee -a ${LOG_FILE}
echo "=========================================" | tee -a ${LOG_FILE}