*IMPORTANT NOTICE*: the UPF format is undergoing a major overhaul and will become a true XML format with a schema. The format described here may become obsolete.

The Unified Pseudopotential Format (UPF), currently at v.2.0.1, is designed to store different kinds of pseudopotentials:

- Norm-conserving (NC) pseudopotentials (PP) in nonlocal form
- As above, in both semilocal (SL) and nonlocal (NL) form
- Ultrasoft (US) PP (aka Vanderbilt)
- Projector Augmented Waves (PAW) datasets
- Additional data for the GIPAW reconstruction of all-electron (AE) charge

into a flexible file format, having a XML-like structure. The PP is stored in numerical form on a radial grid.

## General structure

- The file is formatted, starts with a line
	```
	<UPF version="2.0.1">
	```
	ends with a line
	```
	</UPF>
	```
- The file contains “fields”. A field whose name is “FOO” is delimited by a starting line containing <FOO> and an ending line containing </FOO>, as in the example below:
	```
	<FOO>
	  (content of field FOO)
	</FOO>
	```
	<FOO> and </FOO> are “delimiters” of field FOO.
- A field may have “attributes”, as in the following example:
	```
	<FOO attr="bar">
	  (content of field FOO)
	</FOO>
	```
- A field name is case-sensitive: use all caps in case of doubt!
- A field name can contain only letters and digits
- Spaces are not allowed between the <> and the field name
- A delimiter need not to start at the beginning of a line
- Trailing characters in the line after the > of a delimiter are ignored
- Fields may contain numeric data, character strings, or other fields (“subfields”)
- Blank lines in a field are ignored.
- Comments are introduced by <– and terminated by –> (on the same line)
- Maximum line length is 80 characters
- Data in fields must be readable using fortran free format

## Defined fields

First-level fields defined in v.2.0.1 of the UPF format are:

- PP\_INFO
- PP\_HEADER
- PP\_MESH
- PP\_NLCC (optional)
- PP\_LOCAL
- PP\_NONLOCAL
- PP\_SEMILOCAL (optional, only for norm-conserving)
- PP\_PSWFC (optional)
- PP\_FULL\_WFC (only for PAW)
- PP\_RHOATOM
- PP\_PAW (only for PAW)

PP\_INFO should be the first field, but only for human readers: it is meant to contain info that would allow to generate again the pseudopotential, not data to be read. PP\_HEADER must precede PP\_MESH that must precede all the remaining fields.

Fields that are not defined are ignored.

## Field specifications

The meaning of variables is explained at the end of each field. Do loops, indicated by lines containing only dots (…), can be written with as many numbers per line as desired, within the limit of 80 characters per line.

All quantities are in atomic Rydberg units: e^2=2, m=1/2, hbar=1. Lengths are in Bohr (0.529177 A), energies in Ry (13.6058 eV). Potentials are multiplied by e so they have the units of energy. Beware: some quantities computed on the radial mesh contain a multiplying *r* factor.

### PP\_INFO

The PP\_INFO field may contain any piece of information that is deemed useful to re-generate the pseudopotential. The recommanded structure is the following:

```
<PP_INFO>
  Generated using XXX code v.N
  Author: Jon Doe
  Generation date: 32Oct1976
  Pseudopotential type: SL|NC|1/r|US|PAW
  Element:  Tc
  Functional:  SLA  PW   PBX  PBC
  Suggested minimum cutoff for wavefunctions:  N Ry
  Suggested minimum cutoff for charge density: M Ry
  Non-/scalar-/fully-relativistic pseudopotential
  Local potential generation info (L, rcloc, pseudization)
  Pseudopotential is spin-orbit/contains GIPAW data
  Valence configuration:
  nl, pn, l, occ, Rcut, Rcut US, E pseu
  els(1),  nns(1),  lchi(1),  oc(1),  rcut(1),  rcutus(1),  epseu(1)
  ...
  els(n),  nns(n),  lchi(n),  oc(n),  rcut(n),  rcutus(n),  epseu(n)
  Generation configuration:
     as above, including all states used in generation
  Pseudization used: Martins-Troullier/RRKJ
  <PP_INPUTFILE>
    Copy of the input file used in generation
  </PP_INPUTFILE>
</PP_INFO>
```

---

**PP\_HEADER**

Structure:

```
<PP\_HEADER attr1="value1" ... attrN="valueN"> ... </PP\_HEADER>
     attr             value
  generated       "Generation code"
  author          "Author"
  date            "Generation date"
  comment         "Brief description"
  element         "Chemical Symbol"
  pseudo_type     "NC | SL | 1/r | US | PAW"
  relativistic    "scalar | full | nonrelativistic"
  is_ultrasoft    .F. | .T.
  is_paw          .F. | .T.
  is_coulomb      .F. | .T.
  has_so          .F. | .T.
  has_wfc         .F. | .T.
  has_gipaw       .F. | .T.
  paw\_as\_gipaw    .F. | .T.
  core_correction .F. | .T.
  functional      "dft"
  z_valence        Zval
  total_psenergy   etotps
  wfc_cutoff       ecutwfc
  rho_cutoff       ecutrho
  l_max            lmax
  l\_max\_rho        lmax_rho
  l_local          lloc
  mesh_size        mesh
  number\_of\_wfc    nwfc
  number\_of\_proj   nbeta
</PP_HEADER>
```

---

For cases in which different values of attributes are listed, the first one is the default one. Meaning of variables:

- “NC” = Norm-Conserving PP, fully nonlocal form only
- “SL” = Norm-Conserving PP, nonlocal and semilocal forms available
- “US” = Ultrasoft (Vanderbilt) PP. Implies: is\_uspp=.T.
- “1/r”= Coulomb potential:. Implies: is\_coulomb=.T.
- “PAW”=Projector-Augmented Wave set. Implies: is\_paw=.T
- has\_so: fully relativistic PP with spin-orbit terms
- has\_wfc: contains all-electron orbitaks in field PP\_FULL\_WFC
- has\_gipaw: contains data for GIPAW reconstructions in PP\_GIPAW
- paw\_as\_gipaw=?
- nlcc: non-linear core correction is included
- “dft”: a label identifying the exchange-correlation functional
- Zval: valence charge
- etotps: total pseudo-valence energy of PP
- ecutwfc: suggested cutoff for orbital expansion into plane waves
- ecutrho: suggested cutoff for charge density expansion
- lmax: max angular: momentum component in PP
- lmax\_rho: as above, in atomic charge density (PAW only)
- mesh: number of points in radial grid
- wfc: is the number of atomic (pseudo-)orbitals in section PP\_PSWFC. May not coincide with the number of atomic states used in PP generation
- nbeta: number of Kleinman-Bylander projectors (“beta functions”) included in the PP (field PP\_NONLOCAL)
- - \*

**PP\_MESH**

Structure:

```
<PP_MESH dx="dx" mesh="mesh" xmin="xmin" rmax="rmax" zmesh="Zmesh">
  <PP_R>
     r(1) r(2) ...  r(mesh)
  </PP_R>
  <PP_RAB>
     rab(1) rab(2) ... rab(mesh)
  </PP_RAB>
</PP_MESH>
```

---

Meaning of the variables:

- dx, mesh, xmi, rmax, Zmesh: radial grid arameters
- r (1:mesh): radial grid points (a.u.). Can be one of the following:*  
	r(i) = exmin+i\*dx/Zmesh* or *r(i) = (exmin+i\*dx-1)/Zmesh*, with *r(mesh)=rmax*
- rab(mesh): factor required for discrete integration: *∫ f(r) dr =∑ifi rabi*.

**PP\_NLCC**

```
<PP_NLCC>
  rho\_atc(1) rho\_atc(2) ... rho_atc(mesh)
</PP_NLCC>
```

---

Optional, needed only for PP with core corrections. Meaning of variables:

- rho\_atc(mesh): core charge for nonlinear core correction (true charge, not multiplied by 4π *r2*)

**PP\_SEMILOCAL**

```
<PP_SEMILOCAL>
  <PP_VNL1 L="l1" J="j1">
       V\_1(1) V\_1(2) ... V_1(mesh)
  </PP_VNL1 L="l" J="j">
  <PP_VNL2 L="l2" J="j2">
       V\_2(1) V\_2(2) ... V_2(mesh)
  </PP_VNL2 L="l" J="j">
 ...
</PP_SEMILOCAL>
```

---

Optional, NC PP with semilocal form only. There are nbeta terms, each one containing potential V(1:mesh) for the specified value of l (also j in fully relativistic case only)

**PP\_LOCAL**

```
<PP_LOCAL>
   vloc(1) vloc(2) ... vloc(mesh)
</PP_LOCAL>
```

---

vloc(mesh): local potential (Ry a.u.) sampled on the radial grid

**PP\_NONLOCAL**

Structure:

```
<PP_NONLOCAL>

  <PP_BETA>
    1 lll(1)   "Beta  L"
    kkbeta(1)
    beta(1,1) beta(2,1) ... beta(kkbeta(1),1)
  </PP_BETA>
     ...
  <PP_BETA>
    nbeta lll(nbeta)   "Beta  L"
    kkbeta(nbeta)
    beta(1,nbeta) beta(2,nbeta) ... beta(kkbeta(nbeta),nbeta)
  </PP_BETA>

  <PP_DIJ>
     nd, "Number of nonzero Dij"
do nb=1,nbeta
   do mb=nb,nbeta
      if (abs (dion (nb, mb) ) > 0) then
               nb  mb  dion(nb,mb)    "Q_int"
      end if
   end do
end do
  </PP_DIJ>

  <PP_QIJ>
     nqf   "nqf"
     <PP_RINNER>
        rinner(0) rinner(1) ... rinner(2*lmax)
     </PP_RINNER>
do nb=1,nbeta
   do mb=nb,nbeta
           nb  mb  lll(mb)   "i  j  (l)"
           qqq(nb,mb)    "Q_int"
           qfunc(1, nb,mb) qfunc(2, nb,mb) ... qfunc(mesh, nb,mb)
     <PP_QFCOEF>
           do l=0,2*lmax
              do i=1,nqf
                 qfcoef(i,l,nb,mb)
              end do
           end do
     </PP_QFCOEF>
        end do
     end do
  </PP_QIJ>

</PP_NONLOCAL>
```

---

- lll(i): angular momentum of projector i
- kkbeta(i): number of mesh points for projector i (must be ≤ mesh )
- beta(i): projector *riβ(ri)* (note the factor *r*!)
- dion(i,j): the *Dij* factors of the nonlocal PP: V *NL =* ∑i,j Di,j |βi><βj\_
- nqf: number of expansion coefficients for *qij* (may be zero)
- rinner(i): for *r < rinner(i)* Q functions are pseudized (not read if nqf=0)
- qqq(i,j): *Qij = ∫ qij(r) d3r*
- qfunc: *r2 qij(r)* for *r > rinner(i)*
- qfcoef: expansion coefficients of *r2 qij(r)* for *r < rinner(i)* (not read if nqf=0)
- Note on units: the “beta” and “dion” are defined as in Vanderbilt’s USPP code and should have Bohr^{-1/2} and Ry units, respectively. Since they enter the calculation only as (beta\*D\*beta), some converters may actually produce “dion” in Ry^{-1} units and “beta” in Ry\*Bohr^{-1/2} units instead, as suggested by the Kleinman-Bylader transformation.

**PP\_PSWFC**

```
<PP_PSWFC>
  els(1) lchi(1) oc(1)  "Wavefunction"
  chi(1,1) chi(2,1) ...  chi(mesh,1)
  ..........
  els(natwfc) lchi(natwfc) oc(natwfc)  "Wavefunction"
  chi(1,natwfc) chi(2,natwfc) ... chi(mesh,natwfc)
</PP_PSWFC>
```

---

- chi(mesh,i): *χi(r)*, *i* -th radial atomic (pseudo-)orbital (radial part of the KS equation, multiplied by *r*)
- els(natwf), lchi(natwf), oc(natwf): as in PP\_HEADER

**PP\_RHOATOM**

```
<PP_RHOATOM>
   rho\_at(1) rho\_at(2) ... rho_at(mesh)
</PP_RHOATOM>
```

---

- rho\_at(mesh): radial atomic (pseudo-)charge.This is *4π r2* times the true charge.

## Additional Fields

### PAW

If a PAW dataset is contained in the UPF file then the additional structure <PAW> is present; it contains the fields listed in the following sections.

**PP\_PAW\_FORMAT\_VERSION**

```
<PP\_PAW\_FORMAT_VERSION>
  version number
</PP\_PAW\_FORMAT_VERSION>
```

---

Contains version number, current version is 0.1.

**PP\_AUGMENTATION**

```
<PP_AUGMENTATION>
  Shape of augmentation charge:
  BESSEL | GAUSS | PSQ | ...
  r\_match\_augfun, irc  "augmentation max radius"
  lmax_aug             "augmentation max angular momentum"
  "Augmentation multipoles:"
  nb  = 1,nbeta
    nb1 = 1,nbeta
      l   = 0,lmax_aug
        augmom(nb,nb1,l)
      enddo
    enddo
  enddo
  "Augmentation functions:"
   do l = 0,lmax_aug
       do nb = 1,nbeta
       do nb1 = 1,nbeta
           if (abs(augmom(nb,nb1,l)) > 0)
              "label of the augmentation function"
              augfun(k), k = 1, mesh
           endif
       enddo
       enddo
   enddo
</PP_AUGMENTATION>
```

---

Data required to build augmentation functions:

- BESSEL|GAUSS|PSQ|…: function used to pseudize the augmentation
- r\_match\_augfun: range beyond which all the augmentation functions are zero (a.u.)
- irc: index of the radial grid closer to, and greater than, r\_match\_augfun
- augmom: multipole of augmentation channel (nb,nb1); it is computed as:  
	*mlnb,nb1 =∫ dr r2 rl (χAEnb(r) χAEnb1(r) – χPSnb(r) χPSnb1(r))*
- augfun: the augmentation function, it is stored only if the corresponding augmentation multipole is different from zero.

**PP\_AE\_RHO\_ATC**

```
<PP\_AE\_RHO_ATC>
    aeccharg(k), k = 1,mesh
</PP\_AE\_RHO_ATC>
```

---

All-electron atomic density on the radial grid.

**PP\_AEWFC**

```
<PP_AEWFC>
   do nb = 1,nbeta
      aewfc(k, nb), k = 1,mesh
   end do
</PP_AEWFC>
```

---

All-electron wavefunctions used for the generation of the dataset; there is one wavefunction for each beta projector.

**PP\_PSWFC\_FULL**

```
<PP\_PSWFC\_FULL>
   do nb = 1,nbeta
      pswfc(k, nb), k = 1,mesh
   end do
</PP\_PSWFC\_FULL>
```

---

Pseudo wavefunction used for the generation of the dataset; note that in the PP\_PSWFC field only the occupied wavefunctions are stored while for PAW you need a wavefunction for each projector.

**PP\_AEVLOC**

```
<PP_AEVLOC>
   do nb = 1,nbeta
      aewfc(k, nb), k = 1,mesh
   end do
</PP_AEVLOC>
```

---

All-electron local potential

**PP\_KDIFF**

```
<PP_KDIFF>
  nb  = 1,nbeta
    nb1 = 1,nbeta
      kdiff(nb, nb1)
    enddo
  enddo
<PP_KDIFF>
```

---

Kinetic energy difference between all-electron and pseudo component of each augmentation channel.

**PP\_OCCUP**

```
<PP_OCCUP>
   do nb = 1,nbeta
      occ(nb)
   end do
</PP_OCCUP>
```

---

Occupations of atomic orbitals.

**PP\_GRID\_RECON**

```
<PP\_GRID\_RECON>
   "Minimal info to reconstruct the radial grid:"
   grid%dx,   "  dx"
   grid%xmin, "  xmin"
   grid%rmax, "  rmax"
   grid%zmesh,"  zmesh"

   <PP\_SQRT\_R>
   grid%sqr(k), k=1,mesh
   </PP\_SQRT\_R>
</PP\_GRID\_RECON>
```

---

Addition data necessary to accurately reconstruct the radial grid used for the dataset generation.

**GIPAW**

GIPAW additional data is necessary to reconstruct all-electron charge density using the gipaw.x program included in QE distribution.

**PP\_GIPAW\_FORMAT\_VERSION**

```
<PP\_GIPAW\_FORMAT_VERSION>
  version number
</PP\_GIPAW\_FORMAT_VERSION>
```

---

Contains version number, current version is 0.1.

**GIPAW\_CORE\_ORBITALS**

```
<GIPAW\_CORE\_ORBITALS>
n\_core\_orbitals "number of core orbitals"
    <GIPAW\_CORE\_ORBITAL>
       n, l  "orbital n and l quantum numbers"
       core_orbital(k), k = 1,mesh
    </GIPAW\_CORE\_ORBITAL>
```

Repeated for each core orbital

```
</GIPAW\_CORE\_ORBITALS>
```

---

Core orbitals.

**GIPAW\_LOCAL\_DATA**

```
<GIPAW\_LOCAL\_DATA>
   <GIPAW\_VLOCAL\_AE>
      vlocal_ae(k), k = 1,mesh
   </GIPAW\_VLOCAL\_AE>
   <GIPAW\_VLOCAL\_PS>
      vlocal_ps(k), k = 1,mesh
   </GIPAW\_VLOCAL\_PS>
</GIPAW\_LOCAL\_DATA>
```

---

All electron and pseudo local potentials, sampled on the radial grid.

**GIPAW\_ORBITALS**

```
<GIPAW_ORBITALS>
   <GIPAW\_AE\_ORBITAL>
         el(nb), ll(nb)
         wfs_ae(k,nb), k = 1, mesh
   </GIPAW\_AE\_ORBITAL>
   <GIPAW\_PS\_ORBITAL>
         rcut(nb), rcutus(nb)
         wfs_ae(k,nb), k = 1, mesh
   </GIPAW\_PS\_ORBITAL>
  Repeated for each valence orbital.
</GIPAW_ORBITALS>
```

---

- el: principal quantum number (0,1,2..)
- ll: angular momentum quantum number
- rcut: inner cutof radius (a.u.)
- rcutus: outer cutoff radius (a.u.)
- wfc\_ae: all-electron wavefunction sample on radial grid
- wfc\_ps: pseudo wavefunction sample on radial grid

Last change: July 10, 2018

