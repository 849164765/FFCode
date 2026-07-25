#!/bin/bash
# ============================================================
# 云端 KH 仿真快速诊断 (Bash版, 无需Python)
#
# 用法:
#   cd /path/to/FerroKHinstability/examples/kh2d
#   bash cloud_diagnose.sh [仿真日志文件]
#
# 如果不指定日志文件，会自动查找 *.log, *.out, kh2d_*.o 等
# ============================================================

set -e
set +H  # 禁用 bash 历史扩展避免 ! 字符报错

echo ""
echo "============================================================"
echo "  云端 KH 仿真快速诊断 (Bash版)"
echo "============================================================"
echo ""

# ---- 1. Git 版本检查 ----
echo "[1/3] Git 代码版本检查"
echo "------------------------------------------------------------"
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
# 使用 git rev-parse 查找仓库根目录 (兼容 .git 文件和目录)
PROJ_ROOT=""
SEARCH_DIR="$SCRIPT_DIR"
for i in $(seq 1 8); do
    if GIT_DIR_OUT=$(git -C "$SEARCH_DIR" rev-parse --git-dir 2>/dev/null); then
        PROJ_ROOT=$(git -C "$SEARCH_DIR" rev-parse --show-toplevel 2>/dev/null || echo "$SEARCH_DIR")
        break
    fi
    SEARCH_DIR="$(dirname "$SEARCH_DIR")"
    [ "$SEARCH_DIR" = "/" ] && break
done

if [ -n "$PROJ_ROOT" ]; then
    BRANCH=$(git -C "$PROJ_ROOT" rev-parse --abbrev-ref HEAD 2>/dev/null || echo "unknown")
    HEAD=$(git -C "$PROJ_ROOT" log --oneline -1 2>/dev/null || echo "unknown")
    echo "  项目根目录: $PROJ_ROOT"
    echo "  当前分支:   $BRANCH"
    echo "  当前 HEAD:  $HEAD"

    # 检查是否包含扭曲周期 BC 修复 (最新)
    if git -C "$PROJ_ROOT" log --oneline --all 2>/dev/null | grep -q "a892a8d"; then
        if git -C "$PROJ_ROOT" merge-base --is-ancestor a892a8d HEAD 2>/dev/null; then
            echo "  [OK] 已包含扭曲周期 BC 修复 (a892a8d) - 最新"
        else
            echo "  [FAIL] HEAD 不包含扭曲周期 BC 修复 (a892a8d)"
            echo "    >>> 这是问题根因 (磁场稳定化不足)"
            echo "    >>> 修复: cd $PROJ_ROOT && git pull origin $BRANCH"
            echo ""
            echo "  预期 HEAD:"
            echo "    a892a8d fix: 实现MF扭曲周期边界条件"
        fi
    elif git -C "$PROJ_ROOT" log --oneline --all 2>/dev/null | grep -q "3d84ac1"; then
        if git -C "$PROJ_ROOT" merge-base --is-ancestor 3d84ac1 HEAD 2>/dev/null; then
            echo "  [WARN] 仅包含 Dirichlet BC 修复 (3d84ac1), 缺少扭曲周期 BC (a892a8d)"
            echo "    >>> Dirichlet BC 稳定化仅 9.3%, 扭曲周期 BC 可达 53.1%"
            echo "    >>> 修复: cd $PROJ_ROOT && git pull origin $BRANCH"
        else
            echo "  [FAIL] HEAD 不包含 Dirichlet BC 修复 (3d84ac1)"
            echo "    >>> 这是问题根因"
            echo "    >>> 修复: cd $PROJ_ROOT && git pull origin $BRANCH"
        fi
    else
        echo "  [FAIL] 仓库中未找到任何 BC 修复"
        echo "    >>> 请拉取最新代码: git pull origin $BRANCH"
    fi
else
    echo "  [WARN] 未找到 git 仓库 (可能是 kh0723/ 备份部署, 无 git 信息)"
fi
echo ""

# ---- 1b. kh2d.cpp BC 实现直接检测 (最可靠) ----
echo "[1b/3] kh2d.cpp BC 实现直接检测 (最可靠)"
echo "------------------------------------------------------------"
KH2D_CPP=""
# 查找 kh2d.cpp 文件
for candidate in "$SCRIPT_DIR/kh2d.cpp" "$PROJ_ROOT/examples/kh2d/kh2d.cpp" "$SCRIPT_DIR/../kh2d/kh2d.cpp"; do
    if [ -f "$candidate" ]; then
        KH2D_CPP="$candidate"
        break
    fi
done
# 也检查 kh0723 备份目录
if [ -z "$KH2D_CPP" ] && [ -n "$PROJ_ROOT" ]; then
    for kh0723_cpp in "$PROJ_ROOT/kh0723/examples/kh2d/kh2d.cpp" "$PROJ_ROOT/../kh0723/examples/kh2d/kh2d.cpp"; do
        if [ -f "$kh0723_cpp" ]; then
            KH2D_CPP="$kh0723_cpp"
            echo "  [WARN] 找到 kh0723/ 备份目录的 kh2d.cpp: $kh0723_cpp"
            echo "    >>> 这是问题根因! 云端可能运行旧备份代码"
            break
        fi
    done
fi

if [ -n "$KH2D_CPP" ]; then
    echo "  检测文件: $KH2D_CPP"

    # 检测 Twisted Periodic BC (最新, 正确)
    if grep -q "MF_Per" "$KH2D_CPP" 2>/dev/null && grep -q "Twisted\|扭曲" "$KH2D_CPP" 2>/dev/null; then
        echo "  [OK] 检测到 Twisted Periodic BC (MF_Per + 扭曲校正) - 最新正确版本"
        echo "       H 场应维持 ~100% of H0, 稳定化 73.3%"
    # 检测 Dirichlet BC (次新, 部分修复)
    elif grep -q "psi_ghost=-H0\*x_ghost" "$KH2D_CPP" 2>/dev/null; then
        echo "  [WARN] 检测到 Dirichlet BC (上下壁), 但缺少 Twisted Periodic BC"
        echo "         H 场维持 ~96% H0, 稳定化仅 9.3% (Twisted BC 可达 73.3%)"
        echo "    >>> 修复: git pull origin FerroKHinstability"
    # 检测 buggy Neumann BC (最旧, 有bug)
    elif grep -q "psi_ghost=bPsi.get(id_int)+H0\*vs" "$KH2D_CPP" 2>/dev/null; then
        echo "  [FAIL][FAIL][FAIL] 检测到 buggy Neumann BC (缺 ov 因子)"
        echo "         H 场衰减到 54% H0, 磁场力仅 46% 重力"
        echo "         KH 演化加速 3.84 倍, 20000步即达 t*=6 (应为 76800步)"
        echo ""
        echo "  >>> 这是云端问题的根因! <<<"
        echo "  修复步骤:"
        echo "    1. cd $PROJ_ROOT && git pull origin FerroKHinstability"
        echo "    2. 确认使用 examples/kh2d/ (不是 kh0723/examples/kh2d/)"
        echo "    3. cd examples/kh2d && cp kh2d_backup.ini kh2d.ini"
        echo "    4. make clean && make"
        echo "    5. sbatch dlf.sh (或 mpiexec -n 128 ./kh2d.exe)"
    else
        echo "  [WARN] 无法识别 BC 实现类型, 请手动检查 kh2d.cpp"
        echo "    关键标记:"
        echo "      Twisted BC: grep 'MF_Per' kh2d.cpp  (期望: 多个匹配)"
        echo "      Buggy BC:   grep 'H0\*vs' kh2d.cpp   (期望: 0 匹配)"
    fi

    # 额外检查: kh0723 备份目录存在性
    if [ -n "$PROJ_ROOT" ] && [ -d "$PROJ_ROOT/kh0723" ]; then
        echo ""
        echo "  [FAIL] 检测到 kh0723/ 备份目录存在!"
        echo "    >>> 此目录含 buggy Neumann BC 代码 + dlf.sh 部署脚本"
        echo "    >>> 云端可能误用此目录部署"
        echo "    >>> 建议: rm -rf $PROJ_ROOT/kh0723/  (删除旧备份)"
        echo "             或确保从 examples/kh2d/ 部署"
    fi
else
    echo "  [WARN] 未找到 kh2d.cpp 文件"
fi
echo ""

# ---- 2. kh2d.ini 配置检查 ----
echo "[2/3] kh2d.ini 配置检查"
echo "------------------------------------------------------------"
INI_FILE="$SCRIPT_DIR/kh2d.ini"
if [ -f "$INI_FILE" ]; then
    echo "  配置文件: $INI_FILE"
    for key in Ni Nj Bom TotalStep OutputStep; do
        val=$(grep -E "^\s*${key}\s*=" "$INI_FILE" 2>/dev/null | head -1 | sed 's/.*=\s*//' | sed 's/[;#].*//' | tr -d ' ')
        case $key in
            Ni) expected=512; desc="网格数X (应为512, L=256半域)" ;;
            Nj) expected=512; desc="网格数Y (应为512)" ;;
            Bom) expected=636; desc="磁Bond数" ;;
            TotalStep) expected=80000; desc="总步数" ;;
            OutputStep) expected=2000; desc="输出间隔" ;;
        esac
        if [ "$val" = "$expected" ]; then
            echo "  ✓ $key = $val  ($desc)"
        else
            echo "  ⚠ $key = $val  (期望: $expected)  $desc"
        fi
    done
else
    echo "  ⚠ 未找到 $INI_FILE"
fi
echo ""

# ---- 3. 仿真日志 H 场诊断 ----
echo "[3/3] 仿真日志 H 场诊断"
echo "------------------------------------------------------------"

# 查找日志文件
LOG_FILE="${1:-}"
if [ -z "$LOG_FILE" ]; then
    # 自动查找
    for candidate in "$SCRIPT_DIR"/*.log "$SCRIPT_DIR"/*.out "$SCRIPT_DIR"/kh2d_*.o "$SCRIPT_DIR"/log.txt; do
        if [ -f "$candidate" ]; then
            LOG_FILE="$candidate"
            break
        fi
    done
fi

if [ -n "$LOG_FILE" ] && [ -f "$LOG_FILE" ]; then
    echo "  日志文件: $LOG_FILE"
    echo ""

    # 提取 Hmag 调试行
    # 格式: [DEBUG step=N] Hmag: min=... max=... mean=... range=... (H0=...)
    HLOG=$(grep -c "DEBUG.*Hmag" "$LOG_FILE" 2>/dev/null || echo 0)
    echo "  找到 $HLOG 条 Hmag 调试记录"

    if [ "$HLOG" -gt 0 ]; then
        echo ""
        echo "  前3条记录:"
        grep "DEBUG.*Hmag" "$LOG_FILE" | head -3 | while read line; do
            echo "    $line"
        done
        echo ""
        echo "  最后3条记录:"
        grep "DEBUG.*Hmag" "$LOG_FILE" | tail -3 | while read line; do
            echo "    $line"
        done
        echo ""

        # 提取 H0 和最后的 Hmag mean
        H0=$(grep "DEBUG.*Hmag" "$LOG_FILE" | head -1 | grep -oP 'H0=\K[0-9.e+-]+' || echo "0")
        LAST_HMEAN=$(grep "DEBUG.*Hmag" "$LOG_FILE" | tail -1 | grep -oP 'mean=\K[0-9.e+-]+' | head -1 || echo "0")

        echo "  理论 H0 = $H0"
        echo "  最后 Hmag mean = $LAST_HMEAN"

        # 计算比值 (使用 awk)
        RATIO=$(awk -v h="$LAST_HMEAN" -v h0="$H0" 'BEGIN{if(h0>0) printf "%.4f", h/h0; else print "N/A"}')
        PERCENT=$(awk -v r="$RATIO" 'BEGIN{if(r!="N/A") printf "%.1f", r*100; else print "N/A"}')

        echo "  |H|/H0 (mean) = $RATIO ($PERCENT%)"
        echo ""

        # 诊断
        if [ "$RATIO" != "N/A" ]; then
            IS_GOOD=$(awk -v r="$RATIO" 'BEGIN{if(r>0.90) print 1; else print 0}')
            IS_NEUMANN=$(awk -v r="$RATIO" 'BEGIN{if(r>0.40 && r<0.70) print 1; else print 0}')
            IS_ZERO=$(awk -v r="$RATIO" 'BEGIN{if(r<0.05) print 1; else print 0}')

            if [ "$IS_GOOD" = "1" ]; then
                echo "  [OK][OK][OK] 诊断: H 场维持良好 (>90% of H0)"
                echo "       磁场力应有效稳定 KH 不稳定性"
                echo ""
                echo "  此时 20000 步应对应早期生长阶段 (振幅 ~13-23% L)"
                echo "  不应出现 t*=6 的液桥形态"
                echo "  若仍出现 t*=6, 请检查是否使用扭曲周期 BC (commit a892a8d)"
            elif [ "$IS_NEUMANN" = "1" ]; then
                echo "  [FAIL][FAIL][FAIL] 诊断: Neumann BC (旧版, 有bug)"
                echo "       H 场衰减至 $PERCENT% of H0 (应 >90%)"
                echo "       磁场力不足, KH演化加速约3.84倍"
                echo "       20000步即达t*=6外观 (应为76800步)"
                echo ""
                echo "  >>> 修复方法:"
                echo "      cd $PROJ_ROOT && git pull origin FerroKHinstability"
                echo "      确认 git log -1 显示 commit a892a8d (扭曲周期 BC)"
                echo "      重新编译: make clean && make"
                echo "      重新运行: ./kh2d.exe"
            elif [ "$IS_ZERO" = "1" ]; then
                echo "  [FAIL][FAIL][FAIL] 诊断: H 场几乎为零"
                echo "       MF 格子可能未演化 (Bom<=0 或更早版本代码)"
                echo "       完全无磁场稳定化, KH演化最快"
                echo ""
                echo "  >>> 修复方法:"
                echo "      1. 检查 kh2d.ini 中 Bom=636"
                echo "      2. 拉取最新代码: git pull origin FerroKHinstability"
                echo "      3. 重新编译运行"
            else
                echo "  [WARN] 诊断: H 场异常 ($PERCENT% of H0)"
                echo "       请检查代码版本和配置"
            fi
        fi
    else
        echo "  [WARN] 日志中未找到 Hmag 调试记录"
        echo "    可能是旧版代码(无调试输出)或日志未重定向"
        echo ""
        echo "  建议运行 Python 诊断脚本:"
        echo "    python3 cloud_diagnose.py vtkoutput/vtidata"
    fi
else
    echo "  [WARN] 未找到日志文件"
    echo "    查找位置: $SCRIPT_DIR/*.log, *.out, kh2d_*.o, log.txt"
    echo ""
    echo "  建议:"
    echo "  1. 如果有仿真日志, 运行: bash cloud_diagnose.sh <日志文件>"
    echo "  2. 如果只有VTK输出, 运行: python3 cloud_diagnose.py vtkoutput/vtidata"
fi

echo ""
echo "============================================================"
echo "  诊断完成"
echo "============================================================"
echo ""
