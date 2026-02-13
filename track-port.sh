#!/bin/bash
# Track C++ to C manifold port progress with a loading bar
# Usage: ./track-port.sh [poll_interval_seconds]

INTERVAL=${1:-30}
SRC_C="$(dirname "$0")/src_c"
TARGET_LINES=17000  # approximate C++ source lines to port

RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
CYAN='\033[0;36m'
BOLD='\033[1m'
DIM='\033[2m'
RESET='\033[0m'

bar() {
    local pct=$1 width=40
    local filled=$(( pct * width / 100 ))
    local empty=$(( width - filled ))
    printf "${GREEN}"
    printf '█%.0s' $(seq 1 $filled 2>/dev/null)
    printf "${DIM}"
    printf '░%.0s' $(seq 1 $empty 2>/dev/null)
    printf "${RESET}"
}

prev_commits=0
prev_lines=0
prev_tests=0

while true; do
    clear
    
    # Gather stats
    lines=$(cd "$SRC_C" && cat *.c *.h 2>/dev/null | wc -l | tr -d ' ')
    files=$(cd "$SRC_C" && ls *.c *.h 2>/dev/null | wc -l | tr -d ' ')
    commits=$(cd "$(dirname "$0")" && git --no-pager log --oneline c-port 2>/dev/null | grep "^[a-f0-9]" | grep -c "port:")
    last_commit=$(cd "$(dirname "$0")" && git --no-pager log --oneline c-port -1 2>/dev/null)
    last_time=$(cd "$(dirname "$0")" && git --no-pager log --format="%ar" c-port -1 2>/dev/null)
    
    # Try to get test count from latest binary
    test_count="?"
    if [ -f "$SRC_C/test_manifold" ]; then
        result=$(cd "$SRC_C" && perl -e 'alarm 30; exec("./test_manifold")' 2>&1 | tail -1)
        if echo "$result" | grep -q "passed"; then
            test_count=$(echo "$result" | grep -oE '[0-9]+' | head -1)
        elif echo "$result" | grep -q "FAIL"; then
            test_count=$(echo "$result" | grep -oE '[0-9]+' | head -1)
            test_count="${test_count} (FAILING)"
        fi
    fi
    
    pct=$(( lines * 100 / TARGET_LINES ))
    [ $pct -gt 100 ] && pct=100
    
    # Deltas
    d_commits=""
    d_lines=""
    d_tests=""
    if [ $prev_commits -gt 0 ]; then
        dc=$(( commits - prev_commits ))
        dl=$(( lines - prev_lines ))
        [ $dc -gt 0 ] && d_commits=" ${CYAN}(+${dc})${RESET}"
        [ $dl -gt 0 ] && d_lines=" ${CYAN}(+${dl})${RESET}"
        if [ "$test_count" != "?" ] && [ $prev_tests -gt 0 ]; then
            dt=$(( ${test_count%% *} - prev_tests ))
            [ $dt -gt 0 ] && d_tests=" ${CYAN}(+${dt})${RESET}"
        fi
    fi
    prev_commits=$commits
    prev_lines=$lines
    [[ "$test_count" =~ ^[0-9]+$ ]] && prev_tests=$test_count
    
    echo -e "${BOLD}╔══════════════════════════════════════════════════╗${RESET}"
    echo -e "${BOLD}║     🦎 Manifold C++ → C Port Progress 🦎       ║${RESET}"
    echo -e "${BOLD}╠══════════════════════════════════════════════════╣${RESET}"
    echo ""
    printf "  "
    bar $pct
    printf "  ${BOLD}%d%%${RESET}\n\n" $pct
    echo -e "  ${BOLD}Lines:${RESET}    ${lines} / ${TARGET_LINES}${d_lines}"
    echo -e "  ${BOLD}Files:${RESET}    ${files}"
    echo -e "  ${BOLD}Commits:${RESET}  ${commits}${d_commits}"
    echo -e "  ${BOLD}Tests:${RESET}    ${test_count}${d_tests}"
    echo ""
    echo -e "  ${DIM}Latest:${RESET} ${last_commit}"
    echo -e "  ${DIM}When:${RESET}   ${last_time}"
    echo ""
    echo -e "${BOLD}╚══════════════════════════════════════════════════╝${RESET}"
    echo -e "  ${DIM}Polling every ${INTERVAL}s · Ctrl+C to stop${RESET}"
    
    sleep "$INTERVAL"
done
