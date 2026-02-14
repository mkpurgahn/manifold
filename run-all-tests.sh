#!/bin/bash
# run-all-tests.sh — Sequentially run ralph-loop through all 58 test phases
# Usage: ./run-all-tests.sh [--start-phase N]

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TESTS_DIR="$SCRIPT_DIR/prompts/tests"
RALPH_LOOP="/Users/admin/projects/Axolotl/tools/ralph-loop/ralph-loop.sh"
WORKING_DIR="/Users/admin/projects/manifold"

COMPLETION_PROMISE="Phase complete — test ported and verified"
START_PHASE=1

# Parse args
while [[ $# -gt 0 ]]; do
    case $1 in
        --start-phase)
            START_PHASE="$2"; shift 2 ;;
        -h|--help)
            echo "Usage: ./run-all-tests.sh [--start-phase N]"
            echo "  Runs ralph-loop for each test phase in sequence (1-58)"
            echo "  --start-phase N  Start at phase N (default: 1)"
            exit 0 ;;
        *)
            echo "Unknown arg: $1" >&2; exit 1 ;;
    esac
done

# Collect phase files in order
PHASES=()
for f in "$TESTS_DIR"/[0-9][0-9]-*.md; do
    [[ -f "$f" ]] && PHASES+=("$(basename "$f")")
done

TOTAL=${#PHASES[@]}

echo "═══════════════════════════════════════════════════════════"
echo "🧪 Manifold C Port — Test Phase Sequencer"
echo "═══════════════════════════════════════════════════════════"
echo "Total phases: $TOTAL"
echo "Starting at phase: $START_PHASE"
echo "Working directory: $WORKING_DIR"
echo ""

cd "$WORKING_DIR"

for i in "${!PHASES[@]}"; do
    PHASE_NUM=$((i + 1))

    # Skip phases before start
    if [[ $PHASE_NUM -lt $START_PHASE ]]; then
        echo "⏭️  Skipping Phase $PHASE_NUM (before start phase)"
        continue
    fi

    PHASE_FILE="${PHASES[$i]}"
    PHASE_PATH="$TESTS_DIR/$PHASE_FILE"

    if [[ ! -f "$PHASE_PATH" ]]; then
        echo "❌ Phase file not found: $PHASE_PATH" >&2
        exit 1
    fi

    echo ""
    echo "═══════════════════════════════════════════════════════════"
    echo "🚀 Phase $PHASE_NUM of $TOTAL: $PHASE_FILE"
    echo "═══════════════════════════════════════════════════════════"
    echo ""

    # Run ralph-loop for this phase
    "$RALPH_LOOP" "$PHASE_PATH" \
        --completion-promise "$COMPLETION_PROMISE" \
        --working-dir "$WORKING_DIR"

    EXIT_CODE=$?

    if [[ $EXIT_CODE -ne 0 ]]; then
        echo "❌ Phase $PHASE_NUM failed or was stopped" >&2
        exit $EXIT_CODE
    fi

    echo ""
    echo "✅ Phase $PHASE_NUM complete: $PHASE_FILE"
    echo ""

    # Brief pause between phases
    if [[ $PHASE_NUM -lt $TOTAL ]]; then
        echo "⏳ Pausing 5s before next phase..."
        sleep 5
    fi
done

echo ""
echo "═══════════════════════════════════════════════════════════"
echo "🎉 ALL $TOTAL TEST PHASES COMPLETE!"
echo "═══════════════════════════════════════════════════════════"
