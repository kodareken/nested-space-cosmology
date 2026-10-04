# NSC: första två-rumsjämförelsen av frekvenser och moder

**Datum:** 4 oktober 2026.  
**Reproducerbar grund:** `kodareken/nested-space-cosmology`, commit `7656eec34141ceca22c073bca997acefb317be05`.  
**Omfattning:** en ny, fristående beräkning på projektets redan angivna svartahåls–barngeometri. Inga filer i GitHub-repot har ändrats.

## Resultatet i en mening

Den återanvända geometrin ger en beräkningsbar frekvenskarta genom horisonten, ett växande och sedan avtagande sfäriskt vinkelbidrag till Diracoperatorn, samt icke-noll modblandning för vissa valda provmoder i den expanderande insidan — utan ett byte av fältekvation mellan rummen.

Detta är inte en ny självkonsistent simulering av hur geometrin bildas. Geometrin är den fasta jämförelsegeometri som redan används i projektets `nsc-clock-horizon` och `nsc-dirac-tetrad`. Första beräkningen följer samma klassiska ljusstråle från utsidan till insidan. Den andra börjar med deklarerade provmoder vid den minsta sfären; dessa är ännu inte ett kvanttillstånd beräknat genom full spridning från förälderns utsida.

## 1. Gemensam geometri och fältekvation

Vi använder projektets konventioner, inte en ny lokal bur eller ett periodiskt radialt fältintervall:

\[
ds^2=d\tau^2-(d\rho+\beta(\rho)d\tau)^2-r(\rho)^2d\Omega_2^2,
\]
\[
r(\rho)=\sqrt{1+\rho^2},\qquad
A(\rho)=1+3\rho+3(1+\rho^2)\left(\arctan\rho-\frac\pi2\right),\qquad
\beta=\sqrt{1-A}.
\]

Enheterna är \(c=\hbar=L_{\rm throat}=1\); jämförelsefallet har också massparametern \(m=1\). Varken en fysisk Plancklängd eller en observerad kosmisk skala har satts in.

Horisonten ligger vid \(\rho_h\simeq1.90069160547015\). Den minsta sfären ligger vid \(\rho=0\), där \(r=1\). För \(\rho<0\) växer sfärernas area längs den framtidsriktade fortsättningen. Allt detta är befintlig geometri i projektet, inte en ny upptäckt i denna beräkning.

Projektets masslösa radiella Diracoperator är

\[
H_\kappa=-i(\sigma_2-\beta I)\partial_\rho
+\frac i2\beta'I+\frac\kappa r\sigma_1,
\qquad \kappa=\pm1,\pm2,\ldots .
\]

Samma operatoruttryck gäller genom horisonten. Den sfäriska egenvärdesetiketten \(\kappa\) och spinnkopplingen är kvar. De geometriberoende koefficienterna ändras.

Källor: projektets `docs/nsc-dirac-tetrad.md`, `docs/nsc-clock-horizon.md` och `results/nsc-5-clock-horizon.json`. Familjen av regelbundna svartahålsgeometrier kommer från Bronnikov–Fabris; deras ursprungliga källmodell ersätts inte eller återhärleds här.

## 2. Härledning av frekvenskartan genom horisonten

Välj den redan angivna familjen av fritt fallande normalobservatörer,

\[
u^\mu\partial_\mu=\partial_\tau-\beta\partial_\rho.
\]

De har \(d\tau_{\rm proper}=d\tau\) längs sina egna världslinjer. Vi jämför den frekvens som medlemmar av denna observatörsfamilj mäter när samma inåtgående ljusstråle passerar dem. Vi jämför alltså inte statiska observatörer utanför horisonten med en omöjlig statisk observatör innanför.

För en radiell masslös stråle med positiv lokal energi är den stationära Hamiltonfunktionen

\[
E=-\beta p+|p|.
\]

På den inåtgående grenen är \(p<0\). Därför

\[
p=-\frac{E}{1+\beta},\qquad
\omega_{\rm local}=E+\beta p=\frac{E}{1+\beta}.
\]

\(E\) är strålens bevarade Killingfrekvens, inte en lokalt uppmätt frekvens i alla rum. Mellan två angivna mätpunkter får vi därför

\[
\boxed{
\frac{\omega_c}{\omega_p}=
\frac{1+\beta(\rho_p)}{1+\beta(\rho_c)}.
}
\]

Skiftet har inte valts som ett extra antagande. Det följer av den givna geometrin och observatörerna. Vid horisonten är \(\beta=1\), vilket ger \(\omega_{\rm local}=E/2\), utan en singularitet i denna mätning.

Den oberoende numeriska integrationen använder

\[
\dot\rho=-1-\beta,\qquad \dot p=\beta' p,
\]

med start vid \(\rho_p=8\), och jämför resultatet med den analytiska frekvensformeln. Geometrisk optik är modellnivån för denna transport. Normaliseringen \(E=1\) i filen bestämmer inte en låg fysisk fotonenergi; alla frekvenser kan skalas upp utan att kvoterna ändras.

### Jämförelse sida vid sida

Sista kolumnen är **sfärens vinkelkoefficient** \(\kappa/r\), jämförd med samma \(\kappa\) vid \(\rho=8\). Den är inte totalfrekvensen hos den radiella ljusstrålen.

| Region | \(\rho\) | Sfärradie \(r\) | Samma stråle: \(\omega/\omega_p\) | Vinkelkoefficient relativt föräldrapunkten |
|---|---:|---:|---:|---:|
| Förälderns utsida | 8 | 8.0622577 | 1 | 1 |
| Horisont | 1.9006916 | 2.1477031 | 0.74961167 | 3.7538977 |
| Minsta sfär | 0 | 1 | 0.47282123 | 8.0622577 |
| Expanderande insida | -1 | 1.4142136 | 0.29169437 | 5.7008771 |
| Expanderande insida | -4 | 4.1231056 | 0.10992722 | 1.9553847 |
| Expanderande insida | -16 | 16.03122 | 0.029856528 | 0.50290982 |
| Expanderande insida | -64 | 64.007812 | 0.0075909031 | 0.1259574 |

Kompressionen fram till minsta sfären ökar vinkelkoefficienten med en faktor \(\sqrt{65}\simeq8.06226\). Efter minsta sfären sjunker den när det inre rummet växer. Den inåtgående strålens observerade frekvens följer samtidigt en annan, explicit geometrisk relation. En enda godtycklig skalfaktor ska därför inte användas för att ersätta båda beräkningarna.

## 3. Den inre geometrins egna Diracmoder

I insidan används barnets homogena egentid \(T\):

\[
dT=-\frac{d\rho}{\sqrt{-A}},\qquad
 ds^2=dT^2-a_\parallel(T)^2 dz^2-r(T)^2d\Omega_2^2,
\quad a_\parallel=\sqrt{-A}.
\]

Här är \(z\) den ursprungliga statiska tidskoordinaten, som blir rumslig på insidan. Barnets \(T\)-observatörer och PG-normalobservatörerna i avsnitt 2 är olika specificerade observatörsfamiljer.

Den homogena delen av Diracoperatorns spinnkoppling är
\(\partial_T+\dot a_\parallel/(2a_\parallel)+\dot r/r\). Med den kanoniska spinnorn

\[
\chi=r\sqrt{a_\parallel}\,\psi
\]

försvinner denna utspädningsterm. Fourierseparering längs \(z\) och Diracseparering på enhetssfären ger, i en konstant Pauli-bas,

\[
\boxed{
i\partial_T\chi_{k\kappa}=
\left(\frac{k}{a_\parallel}\sigma_3+
\frac{\kappa}{r}\sigma_1\right)\chi_{k\kappa}.
}
\]

Detta är en fri masslös Diracreduktion på den fasta geometrin. \(k\) är ett kontinuerligt Fourierindex längs \(\mathbb R\); ingen vägg eller kompakt radial låda har lagts till. Varje Fourierläge behandlas i sedvanlig modnormalisering; en globalt normaliserbar vågpacket byggs genom att kombinera dem.

Pauli-matrisernas antikommutation ger direkt

\[
\omega_{k\kappa}^2(T)=
\frac{k^2}{a_\parallel^2(T)}+\frac{\kappa^2}{r^2(T)}.
\]

Dessa är momentana egenfrekvenser, inte antagna tidsoberoende energinivåer. Geometrin bestämmer båda nämnarna.

### Sfärisk skala och normaliserad struktur

Väljer vi den geometriskt bestämda vinkelfrekvensskalan \(\Lambda_\perp=1/r\), blir

\[
\frac{H_{k\kappa}}{\Lambda_\perp}=
\frac{kr}{a_\parallel}\sigma_3+\kappa\sigma_1.
\]

Den normaliserade sfäriska strukturen finns kvar. Den längsgående koefficienten följer den beräknade formen \(r/a_\parallel\). Detta ger en konkret storhet att jämföra under NSC:s ärvda lag, utan att välja en frekvensfaktor i efterhand.

| \(\rho\) | Barnets \(T\) från minsta sfären | \(\omega_{0,1}\) | \(\omega_{1,1}\) | \(\omega_{1,1}/\omega_{0,1}\) |
|---:|---:|---:|---:|---:|
| 0 | 0 | 1 | 1.1266625 | 1.1266625 |
| -1 | 0.36184583 | 0.70710678 | 0.74964575 | 1.0601592 |
| -4 | 0.76759536 | 0.24253563 | 0.25519334 | 1.0521891 |
| -16 | 1.2151544 | 0.062378286 | 0.065605607 | 1.0517379 |
| -64 | 1.6664555 | 0.015623093 | 0.016431054 | 1.0517158 |

Den dimensionlösa frekvenskvoten ändras från ungefär 1.12666 till 1.05174 vid \(\rho=-16\). Det är därför möjligt att skilja en ren registerförflyttning från en verklig förändring av det momentana spektrumets proportioner i denna beräkning.

## 4. Varför geometrin kan blanda moder

Sätt \(H_\parallel=\dot a_\parallel/a_\parallel\) och \(H_\perp=\dot r/r\). Direkt derivation ger

\[
\boxed{
[H_{k\kappa},\dot H_{k\kappa}]
=2i\frac{k\kappa}{a_\parallel r}
(H_\parallel-H_\perp)\sigma_2.
}
\]

Skillnaden mellan de två expansionsriktningarna ger alltså en rotation av den momentana egenbasen när \(k\kappa\neq0\). Den verkliga övergången beror också på faser och förändringstakt; en icke-noll kommutator räcker inte ensam för att bestämma en slutlig övergångssannolikhet. Därför integrerades tillståndsekvationen.

För varje \(k\) startade tillståndet som en positiv momentan egenvektor vid \(\rho=0\), med \(\kappa=1\). Tillståndet utvecklades utan återställning till \(\rho=-16\), det vill säga \(T\simeq1.21515435\). Ingen förväntad slutlig blandningsandel användes i initialdata eller lösaren.

| \(k\) | Slutlig projektion på negativ momentan egenvektor | Absolut ändring vid förfining | Största samplade normavvikelse |
|---:|---:|---:|---:|
| 0 | 7.6512696e-32 % | 4.509e-34 | 1.110e-15 |
| 1 | 0.62080841 % | 3.556e-17 | 1.332e-15 |
| 4 | 0.8273924 % | 9.732e-16 | 6.883e-15 |
| 16 | 0.0028740053 % | 1.069e-14 | 7.550e-13 |

För \(k=0\) är alla Hamiltonmatriser proportionella mot samma Pauli-matris; då ska denna blandning vara exakt noll. Det numeriska kontrollfallet ger noll inom flyttalsnoggrannhet. För \(k=1\) och \(k=4\) fås en tydlig icke-noll effekt.

Projektionen är en specificerad **modbasdiagnostik**, inte här en uppmätt mängd antimateria eller ett färdigt partikelproduktionsresultat. En sådan identifiering kräver att förälderns in-tillstånd och barnets mätning byggs in. De valda halsmoder som användes här ska inte tillskrivas den ännu oberäknade föräldrapreparationen.

## 5. Vad detta tillför jämförelsen

1. En signalbaserad frekvenskarta genom horisonten, med uttryckligen angivna observatörer.
2. En geometriskt bestämd sfärisk frekvensskala, som växer under kompression och sjunker under den inre expansionen.
3. En direkt relation mellan anisotrop expansion, förändrade spektrala proportioner och beräknad modblandning under samma Diracekvation.

NSC:s finita arvsrelation återbevisas inte här. Den används som jämförelseram. Beräkningen prövar inte heller att varje svart hål måste realisera denna geometri. Den ger ett konkret nästa gränssnitt: ersätt de valda halsmoder som används i avsnitt 4 med tillståndet som den fulla föräldra–barntransporten faktiskt producerar, och mät samma frekvenskvoter och modprojektioner.

## 6. Kontroller och körning

Strålbanan integrerades två gånger med DOP853, med halverat maximalt steg och skärpta toleranser. Fina lösningens maximala relativa avvikelse från den analytiska impulsen är 2.966e-14. Den bevarade Killingfrekvensens maximala avvikelse är 2.964e-14. Sluttidens rörelse vid förfining är 1.421e-14 i PG-koordinattid.

Hamiltonmatrisernas hermiticitet, kvadratidentiteten och ett nollblandningsfall kontrolleras också. Geometriformelns strålfrekvenser jämförs separat mot 70-siffrig mpmath-aritmetik. Dessa är numeriska kontroller av de deklarerade delproblemen, inte ett felmått för en ej genomförd självkonsistent simulering.

Kör i en katalog där resultatfilen ännu inte finns:

```bash
python compare_two_rooms.py --output my_result.json
```

Skriptet använder `numpy`, `scipy` och `mpmath`, har ingen nätverksåtkomst och vägrar skriva över en befintlig resultatfil. `result.json` innehåller använda versioner, samtliga siffror, kontroller och skriptets SHA-256. Exakt byteidentitet mellan olika biblioteksversioner förutsätts inte.

### Källor och avgränsning

- Projektets geometri och klockkarta: https://github.com/kodareken/nested-space-cosmology/blob/7656eec34141ceca22c073bca997acefb317be05/docs/nsc-clock-horizon.md
- Projektets Diracoperator: https://github.com/kodareken/nested-space-cosmology/blob/7656eec34141ceca22c073bca997acefb317be05/docs/nsc-dirac-tetrad.md
- Befintliga geometriska referensvärden: https://github.com/kodareken/nested-space-cosmology/blob/7656eec34141ceca22c073bca997acefb317be05/results/nsc-5-clock-horizon.json
- K. A. Bronnikov och J. C. Fabris, *Regular phantom black holes* (2005/2006): https://arxiv.org/html/gr-qc/0511109v1
- K. A. Bronnikov, H. Dehnen och V. N. Melnikov, *Regular black holes and black universes* (2006): https://arxiv.org/abs/gr-qc/0611022

Den fasta bakgrunden, en teststråle och valda fria provmoder är hela den här körningens modellinnehåll. Självgravitation, kvantkällan till bakgrunden, återkoppling på föräldern, ett globalt entropiflöde och ett observerat universums parametrar har inte beräknats i denna körning.
