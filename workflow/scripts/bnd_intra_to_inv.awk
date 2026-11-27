#!/usr/bin/awk -f
# BND intracromosomici con STRANDS ++ o -- -> INV
# Uso: awk -f bnd_intra_to_inv.awk in.vcf > out.vcf
# I record convertiti:
#   SVTYPE=BND -> SVTYPE=INV, ALT -> <INV>, POS = breakpoint minore,
#   END = breakpoint maggiore, SVLEN = END - POS, MATE_ID tolto, flag FROM_BND aggiunto.
# Tutti gli altri record escono invariati.
BEGIN { FS = OFS = "\t"; n_conv = 0; n_bnd = 0 }

/^##/ { print; next }

/^#CHROM/ {
    print "##INFO=<ID=FROM_BND,Number=0,Type=Flag,Description=\"INV ricavata da un BND intracromosomico con STRANDS ++ o --\">"
    print; next
}

{
    if ($8 !~ /(^|;)SVTYPE=BND(;|$)/) { print; next }
    n_bnd++

    # partner dall'ALT: N]chr:pos] , N[chr:pos[ , ]chr:pos]N , [chr:pos[N
    n = split($5, a, /[][]/)
    if (n < 3) { print; next }
    split(a[2], b, ":")
    chr2 = b[1]; pos2 = b[2] + 0

    # STRANDS
    st = ""
    m = split($8, info, ";")
    for (i = 1; i <= m; i++) if (info[i] ~ /^STRANDS=/) st = substr(info[i], 9)

    if (chr2 != $1 || (st != "++" && st != "--")) { print; next }

    p1 = $2 + 0
    lo = (p1 < pos2) ? p1 : pos2
    hi = (p1 < pos2) ? pos2 : p1

    # ricostruisco INFO senza SVTYPE/SVLEN/END/MATE_ID
    newinfo = "SVTYPE=INV"
    for (i = 1; i <= m; i++) {
        if (info[i] ~ /^(SVTYPE|SVLEN|END|MATE_ID)=/) continue
        newinfo = newinfo ";" info[i]
    }
    newinfo = newinfo ";SVLEN=" (hi - lo) ";END=" hi ";FROM_BND"

    $2 = lo; $4 = "N"; $5 = "<INV>"; $8 = newinfo
    n_conv++
    print
}

END { print "BND letti: " n_bnd " | convertiti in INV (intra ++/--): " n_conv > "/dev/stderr" }