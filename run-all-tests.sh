#!/bin/bash
# run-all-tests.sh — Sequentially run ralph-loop through prompt phases
#
# Usage:
#   ./run-all-tests.sh                                    # run all .md in prompts/tests/ (original behavior)
#   ./run-all-tests.sh --dir /path/to/prompts             # run all .md in given directory
#   ./run-all-tests.sh --list prompts.txt                 # run prompts listed in file (one path per line)
#   ./run-all-tests.sh file1.md file2.md file3.md         # run specific prompt files in order
#
# Options:
#   --start-phase N       Start at phase N (1-indexed, default: 1)
#   --dir DIR             Directory of .md prompt files to run (sorted)
#   --list FILE           Text file with one prompt path per line
#   --completion-promise  Override completion promise (auto-extracted from prompt if not set)
#   --working-dir DIR     Working directory for ralph-loop (default: script-detected)
#   --pause SECS          Pause between prompts (default: 5)
#   --dry-run             Print what would run without executing

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DEFAULT_TESTS_DIR="$SCRIPT_DIR/../../prompts/tests"
RALPH_LOOP="/Users/admin/projects/Axolotl/tools/ralph-loop/ralph-loop.sh"
WORKING_DIR="/Users/admin/projects/manifold"

DEFAULT_PROMISE="Phase complete — test ported and verified"
COMPLETION_OVERRIDE=""
START_PHASE=1
PAUSE_SECS=5
DRY_RUN=0
SOURCE_DIR=""
SOURCE_LIST=""

# ── Parse args ──
POSITIONAL=()
while [[ $# -gt 0 ]]; do
    case $1 in
        --start-phase)
            START_PHASE="$2"; shift 2 ;;
        --dir)
            SOURCE_DIR="$2"; shift 2 ;;
        --list)
            SOURCE_LIST="$2"; shift 2 ;;
        --completion-promise)
            COMPLETION_OVERRIDE="$2"; shift 2 ;;
        --working-dir)
            WORKING_DIR="$2"; shift 2 ;;
        --pause)
            PAUSE_SECS="$2"; shift 2 ;;
        --dry-run)
            DRY_RUN=1; shift ;;
        -h|--help)
            head -19 "$0" | tail -18
            exit 0 ;;
        *)
            POSITIONAL+=("$1"); shift ;;
    esac
done

# ── Extract completion promise from a prompt file ──
extract_promise() {
    local file="$1"
    local promise
    promise=$(grep -oE 'Say \*\*"[^"]+"' "$file" 2>/dev/null | head -1 | sed 's/Say \*\*"//;s/"$//')
    if [[ -z "$promise" ]]; then
        promise=$(grep -A2 "Completion Promise" "$file" 2>/dev/null | grep -oE '"[^"]+"' | head -1 | tr -d '"')
    fi
    echo "$promise"
}

# ── Collect prompt files ──
PHASES=()

if [[ ${#POSITIONAL[@]} -gt 0 ]]; then
    # Explicit file arguments
    for f in "${POSITIONAL[@]}"; do
        if [[ -f "$f" ]]; then
            PHASES+=("$f")
        else
            echo "Warning: skipping '$f' (not found)" >&2
        fi
    done
elif [[ -n "$SOURCE_DIR" ]]; then
    # Directory mode
    for f in "$SOURCE_DIR"/*.md; do
        [[ -f "$f" ]] && PHASES+=("$f")
    done
elif [[ -n "$SOURCE_LIST" ]]; then
    # List file mode
    while IFS= read -r line; do
        line="${line%%#*}"
        line="$(echo "$line" | xargs)"
        [[ -n "$line" && -f "$line" ]] && PHASES+=("$line")
    done < "$SOURCE_LIST"
else
    # Default: original behavior — prompts/tests/ directory
    TESTS_DIR="$DEFAULT_TESTS_DIR"
    [[ ! -d "$TESTS_DIR" ]] && TESTS_DIR="$SCRIPT_DIR/prompts/tests"
    [[ ! -d "$TESTS_DIR" ]] && { echo "Error: no prompt source specified and default tests dir not found." >&2; exit 1; }
    for f in "$TESTS_DIR"/[0-9][0-9]-*.md; do
        [[ -f "$f" ]] && PHASES+=("$(basename "$f")")
    done
    # Convert basenames back to full paths
    FULL_PHASES=()
    for p in "${PHASES[@]}"; do
        FULL_PHASES+=("$TESTS_DIR/$p")
    done
    PHASES=("${FULL_PHASES[@]}")
fi

if [[ ${#PHASES[@]} -eq 0 ]]; then
    echo "Error: no prompt files found." >&2
    exit 1
fi

TOTAL=${#PHASES[@]}

# ── Banner ──
echo "═══════════════════════════════════════════════════════════"
echo "🔄 Prompt Phase Sequencer"
echo "═══════════════════════════════════════════════════════════"
echo "Total phases: $TOTAL"
echo "Starting at phase: $START_PHASE"
echo "Working directory: $WORKING_DIR"
[[ -n "$COMPLETION_OVERRIDE" ]] && echo "Promise override: $COMPLETION_OVERRIDE"
echo ""

# ── List phases ──
for i in "${!PHASES[@]}"; do
    NUM=$((i + 1))
    NAME="$(basename "${PHASES[$i]}")"
    if [[ $NUM -lt $START_PHASE ]]; then
        echo "  ⏭️  $NUM. $NAME (skip)"
    else
        P="$COMPLETION_OVERRIDE"
        [[ -z "$P" ]] && P="$(extract_promise "${PHASES[$i]}")"
        [[ -z "$P" ]] && P="$DEFAULT_PROMISE"
        echo "  📋 $NUM. $NAME → \"$P\""
    fi
done
echo ""

[[ $DRY_RUN -eq 1 ]] && { echo "(dry run — nothing executed)"; exit 0; }

# ── Execute ──
cd "$WORKING_DIR"

PASSED=0
FAILED=0

for i in "${!PHASES[@]}"; do
    PHASE_NUM=$((i + 1))
    PHASE_PATH="${PHASES[$i]}"
    PHASE_NAME="$(basename "$PHASE_PATH")"

    if [[ $PHASE_NUM -lt $START_PHASE ]]; then
        echo "⏭️  Skipping Phase $PHASE_NUM (before start phase)"
        continue
    fi

    if [[ ! -f "$PHASE_PATH" ]]; then
        echo "❌ Phase file not found: $PHASE_PATH" >&2
        FAILED=$((FAILED + 1))
        continue
    fi

    # Resolve completion promise
    PROMISE="$COMPLETION_OVERRIDE"
    if [[ -z "$PROMISE" ]]; then
        PROMISE="$(extract_promise "$PHASE_PATH")"
    fi
    if [[ -z "$PROMISE" ]]; then
        PROMISE="$DEFAULT_PROMISE"
    fi

    echo ""
    echo "═══════════════════════════════════════════════════════════"
    echo "🚀 Phase $PHASE_NUM of $TOTAL: $PHASE_NAME"
    echo "   Promise: \"$PROMISE\""
    echo "═══════════════════════════════════════════════════════════"
    echo ""

    "$RALPH_LOOP" "$PHASE_PATH" \
        --completion-promise "$PROMISE" \
        --working-dir "$WORKING_DIR"

    EXIT_CODE=$?

    if [[ $EXIT_CODE -ne 0 ]]; then
        echo "❌ Phase $PHASE_NUM failed: $PHASE_NAME (exit $EXIT_CODE)" >&2
        FAILED=$((FAILED + 1))
        continue
    fi

    PASSED=$((PASSED + 1))
    echo ""
    echo "✅ Phase $PHASE_NUM complete: $PHASE_NAME"
    echo ""

    # Brief pause between phases
    if [[ $PHASE_NUM -lt $TOTAL ]]; then
        echo "⏳ Pausing ${PAUSE_SECS}s before next phase..."
        sleep "$PAUSE_SECS"
    fi
done

echo ""
echo "═══════════════════════════════════════════════════════════"
echo "📊 Results: $PASSED passed, $FAILED failed out of $TOTAL"
if [[ $FAILED -eq 0 ]]; then
    echo "🎉 ALL $TOTAL PHASES COMPLETE!"
fi
echo "═══════════════════════════════════════════════════════════"

[[ $FAILED -gt 0 ]] && exit 1
exit 0
