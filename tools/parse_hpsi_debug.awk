BEGIN {
    printf "%-4s %5s  %14s %14s %14s %14s %14s %10s\n",
        "ikpt", "iband", "<T>", "<V_loc>", "<V_NL>", "E_file", "E_calc", "DE"
    ik = -1
    T=""; V=""; N=""; E=""
}

# block header
$1 == "[HPSI_DEBUG]" && /ikpt=.*iband/ {
    sub(/ikpt=/,        "", $2); ik = int($2)
    sub(/iband\(WG\)=/, "", $3); ib = int($3)
    T=""; V=""; N=""; E=""
}

# expectation values
$1 == "[HPSI_DEBUG]" && index($0, "|T|")   { T = $4 }
$1 == "[HPSI_DEBUG]" && index($0, "V_loc") { V = $4 }
$1 == "[HPSI_DEBUG]" && index($0, "V_NL")  { N = $4 }
$1 == "[HPSI_DEBUG]" &&  $2 == "E_total"   { E = $4 }

# test data line:  ikpt  iband  E_file  E_calc  ...
$1 ~ /^[0-9]+$/ && $2 ~ /^[0-9]+$/ && $0 !~ /^\[/ {
    if ($1+0 == ik && $2+0 == ib && T != "" && E != "") {
        ef = $3 + 0
        printf "%-4d %5d  %14.8f %14.8f %14.8f %14.8f %14.8f %10.6f\n",
            ik+1, ib+1, T+0, V+0, N+0, ef, E+0, ef - E
    }
}
