#!/usr/bin/env bash

set -u

green=''
red=''
yellow=''
blue=''
reset=''

pass_count=0
fail_count=0

print_section() {
	printf '\n%s== %s ==%s\n' "$blue" "$1" "$reset"
}

ok() {
	printf '%s[OK]%s %s\n' "$green" "$reset" "$1"
	pass_count=$((pass_count + 1))
}

bad() {
	printf '%s[FAIL]%s %s\n' "$red" "$reset" "$1"
	fail_count=$((fail_count + 1))
}

have_cmd() {
	command -v "$1" >/dev/null 2>&1
}

show_system_info() {
	print_section "System Info"

	if [ -r /etc/os-release ]; then
		. /etc/os-release
		printf 'OS: %s %s\n' "${PRETTY_NAME:-Unknown}" "${VERSION:-}"
	else
		printf 'OS: %s\n' "$(uname -srm)"
	fi

	printf 'Kernel: %s\n' "$(uname -r)"
	printf 'Host: %s\n' "$(hostname 2>/dev/null || echo unknown)"
	printf 'CPU: %s\n' "$(lscpu 2>/dev/null | awk -F: '/Model name/ {gsub(/^ +/, "", $2); print $2; exit}')"
	printf 'CPU Cores: %s\n' "$(nproc 2>/dev/null || printf 'unknown')"
	printf 'RAM: %s\n' "$(free -h 2>/dev/null | awk '/Mem:/ {print $2 " total, " $3 " used, " $4 " free"}')"
	printf 'Storage: %s\n' "$(df -h / 2>/dev/null | awk 'NR==2 {print $2 " total, " $3 " used, " $4 " avail on " $6}')"

	if have_cmd nvidia-smi; then
		printf 'GPU: %s\n' "$(nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null | paste -sd ', ' -)"
	elif have_cmd lspci; then
		gpu_list=$(lspci 2>/dev/null | awk -F': ' '/VGA compatible controller|3D controller|Display controller/ {print $2}')
		if [ -n "${gpu_list:-}" ]; then
			printf 'GPU: %s\n' "$gpu_list"
		else
			printf 'GPU: none detected\n'
		fi
	else
		printf 'GPU: unable to detect (no nvidia-smi/lspci)\n'
	fi
}

check_cmd_version() {
	local name=$1
	local version_cmd=$2
	local min_version=${3:-}

	if ! have_cmd "$name"; then
		bad "$name not found"
		return
	fi

	if [ -n "$min_version" ] && have_cmd python3; then
		current_version="$($version_cmd 2>/dev/null | head -n1)"
		printf '%s[OK]%s %s found (%s)\n' "$green" "$reset" "$name" "${current_version:-version unavailable}"
	else
		printf '%s[OK]%s %s found\n' "$green" "$reset" "$name"
	fi
	pass_count=$((pass_count + 1))
}

check_python_module() {
	local module=$1
	local label=${2:-$1}

	if ! have_cmd python3; then
		bad "python3 not found; cannot check Python module $label"
		return
	fi

	if python3 - <<PY >/dev/null 2>&1
import importlib.util
raise SystemExit(0 if importlib.util.find_spec('${module}') else 1)
PY
	then
		ok "Python module $label found"
	else
		bad "Python module $label missing"
	fi
}

check_requirements() {
	print_section "WQ Build Requirements"

	if have_cmd cmake; then
		cmake_version=$(cmake --version 2>/dev/null | head -n1)
		ok "cmake found (${cmake_version:-version unavailable})"
	else
		bad "cmake not found"
	fi

	if have_cmd gfortran; then
		ok "gfortran found ($({ gfortran --version 2>/dev/null | head -n1; } | sed 's/^/version /'))"
	else
		bad "gfortran not found"
	fi

	if have_cmd mpirun; then
		ok "mpirun found ($({ mpirun --version 2>/dev/null | head -n1; } | sed 's/^/version /'))"
	else
		bad "mpirun not found"
	fi

	if have_cmd mpifort; then
		ok "mpifort found"
	elif have_cmd mpif90; then
		ok "mpif90 found (usable MPI Fortran compiler)"
	else
		bad "MPI Fortran compiler not found (mpifort/mpif90)"
	fi

	if have_cmd python3; then
		ok "python3 found ($(python3 --version 2>&1))"
	else
		bad "python3 not found"
	fi

	check_python_module numpy
	check_python_module matplotlib

	if have_cmd doxygen; then
		ok "doxygen found (optional docs target)"
	else
		printf '%s[INFO]%s doxygen not found (optional)\n' "$yellow" "$reset"
	fi
}

main() {
	printf '%sWQ Requirement Check%s\n' "$blue" "$reset"
	show_system_info
	check_requirements

	print_section "Summary"
	printf 'Passed: %d\n' "$pass_count"
	printf 'Failed: %d\n' "$fail_count"

	if [ "$fail_count" -eq 0 ]; then
		printf '%sAll required build checks passed.%s\n' "$green" "$reset"
		exit 0
	fi

	printf '%sSome required build checks are missing.%s\n' "$red" "$reset"
	exit 1
}

main "$@"
