#!/usr/bin/env bash
# VarForge Task Loop
# Executes tasks sequentially, respecting dependencies.
# Each task is implemented by spawning a Claude subagent with the task context.
#
# Usage:
#   ./scripts/task-loop.sh                    # Run all tasks from beginning
#   ./scripts/task-loop.sh --start 5          # Start from task 5
#   ./scripts/task-loop.sh --only 3           # Run only task 3
#   ./scripts/task-loop.sh --phase 1          # Run all Phase 1 tasks
#   ./scripts/task-loop.sh --check            # Check which tasks are done

set -euo pipefail

TASK_DIR="docs/plan/tasks"
STATUS_FILE="docs/plan/task-status.json"

# Phase definitions
PHASE_1=(1 2 3 4)
PHASE_2=(5 6 7)
PHASE_3=(8 9 10)
PHASE_4=(11 12 13 14)
PHASE_5=(15 16 17 18)
PHASE_6=(19 20 21)

# Dependencies (task -> prerequisite tasks)
declare -A DEPS
DEPS[1]=""
DEPS[2]=""
DEPS[3]=""
DEPS[4]=""
DEPS[5]="1"
DEPS[6]="1 2 3 4 5"
DEPS[7]="5 6"
DEPS[8]="5"
DEPS[9]="5"
DEPS[10]="5"
DEPS[11]="5"
DEPS[12]="5"
DEPS[13]="6"
DEPS[14]="6"
DEPS[15]="6"
DEPS[16]="13"
DEPS[17]="12"
DEPS[18]="5"
DEPS[19]="6"
DEPS[20]="7"
DEPS[21]=""

# Initialize status file if missing
if [ ! -f "$STATUS_FILE" ]; then
    echo '{}' > "$STATUS_FILE"
fi

get_status() {
    local task_num=$1
    python3 -c "
import json
with open('$STATUS_FILE') as f:
    s = json.load(f)
print(s.get('task_$task_num', 'pending'))
" 2>/dev/null || echo "pending"
}

set_status() {
    local task_num=$1
    local status=$2
    python3 -c "
import json
with open('$STATUS_FILE') as f:
    s = json.load(f)
s['task_$task_num'] = '$status'
with open('$STATUS_FILE', 'w') as f:
    json.dump(s, f, indent=2)
"
}

check_deps() {
    local task_num=$1
    local deps="${DEPS[$task_num]}"
    if [ -z "$deps" ]; then
        return 0
    fi
    for dep in $deps; do
        local dep_status=$(get_status $dep)
        if [ "$dep_status" != "done" ]; then
            echo "Task $task_num blocked: dependency task $dep is $dep_status"
            return 1
        fi
    done
    return 0
}

check_all() {
    echo "=== VarForge Task Status ==="
    for i in $(seq 1 21); do
        local padded=$(printf "%02d" $i)
        local status=$(get_status $i)
        local task_file="$TASK_DIR/task-${padded}-*.md"
        local task_name=$(basename $task_file .md 2>/dev/null | sed 's/task-[0-9]*-//' || echo "unknown")
        printf "  Task %2d [%-8s] %s\n" "$i" "$status" "$task_name"
    done
}

run_task() {
    local task_num=$1
    local padded=$(printf "%02d" $task_num)
    local task_file=$(ls $TASK_DIR/task-${padded}-*.md 2>/dev/null | head -1)

    if [ -z "$task_file" ]; then
        echo "ERROR: No task file found for task $task_num"
        return 1
    fi

    local status=$(get_status $task_num)
    if [ "$status" = "done" ]; then
        echo "Task $task_num already done, skipping."
        return 0
    fi

    if ! check_deps $task_num; then
        echo "Skipping task $task_num due to unmet dependencies."
        return 1
    fi

    echo ""
    echo "============================================"
    echo "  Starting Task $task_num: $(basename $task_file .md)"
    echo "============================================"
    echo ""

    set_status $task_num "in_progress"

    # The task file contains all context needed for a subagent
    echo "Task file: $task_file"
    echo "Read this file and execute the task."
    echo ""
    echo "To run with Claude Code:"
    echo "  claude --model sonnet \"Read $task_file and implement everything it describes. Write tests. Run cargo test and cargo clippy to verify. Commit when done.\""
    echo ""

    return 0
}

# Parse arguments
START=1
ONLY=""
PHASE=""
CHECK=false

while [[ $# -gt 0 ]]; do
    case $1 in
        --start) START=$2; shift 2 ;;
        --only) ONLY=$2; shift 2 ;;
        --phase) PHASE=$2; shift 2 ;;
        --check) CHECK=true; shift ;;
        *) echo "Unknown option: $1"; exit 1 ;;
    esac
done

if $CHECK; then
    check_all
    exit 0
fi

if [ -n "$ONLY" ]; then
    run_task $ONLY
    exit $?
fi

if [ -n "$PHASE" ]; then
    phase_var="PHASE_${PHASE}[@]"
    for task_num in "${!phase_var}"; do
        run_task $task_num || true
    done
    exit 0
fi

# Run all tasks from START
for task_num in $(seq $START 21); do
    run_task $task_num || true
done

echo ""
echo "=== Loop complete ==="
check_all
