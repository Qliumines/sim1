#!/usr/bin/env bash

ASSIGNNAME="criticalpath"
EXPERS="1201 161009 171717 133662 leemis79" 
METRICS="cp"
GOLD_SEED=182330
GOLD_RUNS=10000
EXAMPLE_SEEDS="997431 64138225 80755 268520 70182625"
EXAMPLE_N=5
EXAMPLE_RUNS=300
GRADE_SEED=142177
GRADE_RUNS=1000

ascii_sort()
{
	sort
}

filter_prob_one()
{
	sed -e 's/^[[:space:]]*:[^[:space:]]*:[[:space:]]*//' |\
		sed -e 's/\.0*//' | \
		sed -e 's/[eE][-+]\{0,1\}0\{1,\}[[:space:]]*//' |
	grep '^1$' | tr -d '\n'
}

experiment()
{
	local theSIM="${1}"
	local NET="${2}"
	local RUNS="${3}"
	local SEED="${4}"		
	local missingdata=0
	local residfmt="${5:-__residual-%s-%s-%d.%s}"
	local residrand
	local residoutput
	local output
	local gp
	local missingdata=0

	local netbn="${NET##*/}"
	local exper="${netbn#san-}"
	local exper="${exper%.net}"

	residrand=`printf ${5:-${residfmt}} "$exper" random 0 dat 2>/dev/null`
	residoutput=`printf ${5:-${residfmt}} "$exper" output 0 log 2>/dev/null`
	output=_cp-${exper}-${SEED}-${RUNS}.dat
	eszero "${theSIM}" <(Random ${SEED}| tee "${residrand}") "${RUNS}" "${NET}" |\
		tee "${residoutput}" |\
		"${SIMGRADING}/output-line" |\
		ascii_sort >"${output}"

	test -f "${output}" -a ! -s "${output}" && missingdata=1
	return ${missingdata}
}

makeplotdata()
{
	local gold="${1}"
	local output="${2}"
	paste <(awk '{print $2;}' "${gold}") \
		  <(awk '{print $2;}' "${output}") > "_${output}"
}


###
# SIMS and usage should be defined before sourcing sim-lib.sh
###
SIMS=criticalpath
usage()
{
	cat <<EoT
ASSIGNMENT SPECIFIC OPTIONS 

  $ /assignment/specific/path/grader.sh . [RUNS [SEED]]

Where

  . is the SIM location (required first parameter)

  RUNS is the number of RUNS to use for SIM execution and plot generation

  SEED is a positive integer for a SEED to use for random file inputs

To retain residual data files to (maybe) assist in debugging, export

  $ export GRADER_SAVE_RESIDUALS=${SIMS}

EoT
}

# if SIMGRADING is unknown
test -z "${SIMGRADING}" -a -r ~khellman/SIMGRADING/sim-lib.sh && ~khellman/SIMGRADING/setup.sh ~khellman/SIMGRADING
if ! test -r "${SIMGRADING}/sim-lib.sh" ; then
	cat >&2 <<EoT
ERROR:  SIMGRADING is not in your environment or SIMGRADING/sim-lib.sh cannot
be found.  Have you followed the grader tarball setup.sh instructions on the
assignment's Wiki page?
EoT
	exit 1
fi

test -n "${SIMGRADING}" && source "${SIMGRADING}/sim-lib.sh"

set -e

RUNS=${1:-${GRADE_RUNS}}
SEED=${2:-${GRADE_SEED}}


test_nonexist_tracefile "${simloc}/SIM" MISSINGRANDOM 1 "${graderloc}/simple-san.net"

test_truncated_tracefile "${simloc}/SIM" TRUNCRANDOM 100 "${graderloc}/san-leemis79.net"

# check that edge weights are assigned once per replication, *not* each time
# an edge is traversed during path sum accumulation. 
# replication_allocation.net has an a1/2 edge traversed for each of three
# possible paths, we use an N of 1 and provide only 
grader_msg <<EoT
Testing edge weight allocation algorithm, if the result is the same as a
truncated trace file, then the WRONG algorithm is being used for edge weight
allocation.  Should be  exit status 0  OUTPUT grep 1.
This fails when edge weights are determined on each edge traversal within a
single replication, instead of once per replication.
EoT
export test_truncated_tracefile_evaltest=_test_exit_eq_0_with_output
export _test_output_exitstatus_showexpected=0
test_truncated_tracefile "${simloc}/SIM" TRUNCRANDOM_7  1 "${graderloc}/replication_allocation.net"
unset test_truncated_tracefile_evaltest
unset _test_output_exitstatus_showexpected


# check that multiple critical paths can be found
# "improbable errors" are STILL errors!
mcotmp=`grader_mktemp multiple_critical output`
grader_msg <<EoT
Testing that multiple critical paths can be found with
  "${simloc}/SIM" <(yes 0.5) 1000 "${graderloc}/multiple_critical.net"
EoT
eszero "${simloc}/SIM" <(yes 0.5) 1000 "${graderloc}/multiple_critical.net" |\
	"${SIMGRADING}/output-line" |\
	ascii_sort >"${mcotmp}"
if test x`filter_prob_one <"${mcotmp}"` != "x11"; then
	grader_msg <<EoT
'${simloc}/SIM' fails the automated multiple critical paths test...
EoT
	cat "${mcotmp}"
	grader_msg <<EoT
Confirm there are TWO critical paths with probability 1.0000 in
the output above (they should be  a1/3,a3/5  and  a1/4,a4/5).
EoT
	grader_keystroke
else 
	grader_msg <<EoT
... multiple critical paths passed by simple output comparison.
EoT
fi



# finally, critical path results
declare -a nets
missingdata=0
for NET in `(cd "${graderloc}" && ls -1 san-*.net)` ; do 
	
	netbn="${NET##*/}"
	exper="${netbn#san-}"
	exper="${exper%.net}"
	nets[${#nets[@]}]="${exper}"

	gold="${graderloc}/cp-${exper}-${GOLD_SEED}-${GOLD_RUNS}.dat"

	grader_msg <<EoT
Running SIM on ${NET} ...
  "${simloc}/SIM" <(Random ${SEED}) "${RUNS}" "${NET}"
EoT
	# will leave behind _cp-exper-RUNS-SEED.dat
	output=_cp-${exper}-${SEED}-${RUNS}.dat
	missingdata=0
	experiment "${simloc}/SIM" "${graderloc}/${NET}" ${RUNS} ${SEED} || missingdata=$?
	# check results, generate pdf graph for this SAN
	if ! diff <(awk '{print $1}' "${gold}") <(awk '{print $1}' "${output}"); then
		grader_msg <<EoT
ERROR: number of paths or path description strings are wrong for ${netbn}.

The output above is the result of diff(1).  Lines beginning with < are paths your
SIM is not reporting, lines beginning with > are paths your SIM is reporting, but
that don't actually exist in the network.
EoT
		grader_keystroke
	fi
	# confirm significant digit and scientific notation requirement
	cat "${output}" | while read path prob ; do 
		if ! echo "${prob}" | test_significant_digits 6 >/dev/null ; then
		grader_msg <<EoT
ERROR: not all probabilities for ${netbn} have proper number of significant digits.
  ${path} ${prob}
EoT
			grader_keystroke 
			break
		fi
	done

	makeplotdata "${gold}" "${output}"
	gplotpdf cp ${exper} ${RUNS} ${SEED} ${RUNS} 

	grader_cleanup_experiment ${missingdata} "${METRICS}" "${exper}" "${RUNS}" "${SEED}"
    grader_save_residuals ${exper}
done

# these are the experiments that were run
EXPERS="${nets[*]}"
grader_signoff ${RUNS} ${SEED}

