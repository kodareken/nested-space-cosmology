# NSC: omkopplad tillståndslag → stängd lokal incoming-gate → citerbar publikation

## 1. Bekräftat slutmål, avgränsning och verifierat utgångsläge

**Douglas har bekräftat sprintens avgränsning: lokal incoming-gate. Ingen fråga om lokal kontra global gate återstår.**

Sprinten avslutas när den omkopplade tillståndslagen har **ett nytt, ärligt lokalt gate-resultat**, med reproducerbara register, en läsbar PDF och en offentlig GitHub-version som innehåller resultatets hela beroendekedja.

Gaten gäller båda incoming-constrainterna på samma sammanhängande intervall \(I\) med positiv längd, under

\[
C_\Sigma[g]=U_gC_{\rm up}U_g^\dagger,
\]

med oförändrad upstream-källa och full retarderad \(\delta C_\Sigma[g]\). Ett stängt resultat är antingen lokal EXISTENCE med residualer och felbudget inom tolerans, eller scoped NON-EXISTENCE med en nödvändig relation som utesluter den deklarerade lokala historikklassen.

**Global parent/child-matching, extended stationarity, metric timestep och observationskedjan är efterföljande forskning. De är inte sprintkrav och får inte införas som dolda slutkrav.**

**Abstract, README och release ska uttryckligen kalla resultatet lokalt**, ange intervallet och historikklassen samt hålla den lokala slutsatsen skild från global kosmologisk slutning.

Den genomförda granskningen har följt kod, resultatregister, hashbindningar, Git-historik och publiceringsverktyg. Två Grok-utredningar och en riktad uppföljning har använts. Granskningen startade inga vetenskapliga generatorer och ändrade inga spårade filer.

| Del | Verifierat utgångsläge från granskningen |
|---|---|
| Laboratoriet | `codex/nsc-closure-verification`, tip `5cf59e0`. |
| Omkopplingen | `d5645bd` äger evolverat incoming-tillstånd; `19505e5` kopplar in matter och retarderad derivata i båda constrainterna. |
| Derivatakontrollen | Residual `1.3523e-9` mot kontrolltolerans `3e-8`. Detta är verifierad diskret koppling, ännu ingen fysisk rot. |
| Numerisk värdenoggrannhet | Exakt källfas minskar referensdriften cirka 590 gånger. Återstående matter-drift är cirka `8.9e-6` och `7.5e-6`; felgrinden är fortfarande öppen. |
| Källnoggrannhet | V5 certifierar täckta bakgrundsregioner. Återstående low/subgap och den ändrade historiens bidrag saknar fullständig felbudget. |
| Offentlig version | Separat checkout på `/Users/admin/Documents/nested-space-cosmology`, samma tip `37a3629` som GitHub vid granskningen. Där finns den hashverifierade PDF:n `0.26.0`, 43 sidor. |

Kopplingen dokumenteras i [kodkartan](/Users/admin/Documents/BlackHoles-Infinity/docs/active-code-map.md:3) och [constraintägaren](/Users/admin/Documents/BlackHoles-Infinity/docs/nsc-evolved-incoming-constraints.md:1).

Följande är **låst återbruk — öppnas inte på nytt som reparerbara fel**:

- Homogen, frekvensdiagonal KS-klass utan interface: NON-EXISTENCE med dess registrerade shift-residual cirka `−0.0135927`.
- Slät kompakt parent-historia med exakt fryst incoming `C0`: NON-EXISTENCE. `P_K < −0.022178625` är en strikt gräns, inte det exakta värdet på strömmen.
- Compensator-, första ordningens och jet-rigiditetsresultaten i den frysta klassen.
- De certifierade geometriska koefficienterna, principalmatrisen och den villkorliga lokala ytexistensen. De används som matematiska byggdelar.
- Certifierade bakgrundsregioner genom oändligheten och group32-fönstret `[24,32]`.
- Den deklarerade actionbokföringen: Einstein–GHY, lokala seam-termer och hela Gaussian CTP räknas en gång. `Γ_rest` blir ingen ny fri kraft.

De frysta `C0`-resultaten är regressioner. De uppfyller inte kravet på ett nytt gate-resultat under \(C_\Sigma[g]\).

## 2. Genomförande med leverans efter varje fas

### Fas A — Gör den omkopplade utvärderingen tillförlitlig i värde

Återanvänd samma Dirac-generator, källkolumner, retarderade variationsoperator och restriktion till incoming-ytan. Behåll

\[
C_\Sigma[g]=F[g]C_{\rm src}F[g]^\dagger,\qquad
\delta C_{\rm src}=0.
\]

Inför en sammanhållen produktionsväg som returnerar **\(F,\delta F,F_z,\delta F_z\)** från samma evolution. Källans harmoniska faser evolveras i det utökade systemet. Axialderivatorna tas från den faktiska PDE:n vid restriktionen; källenergin används fortsatt som etikett.

Den kvarlämnade boundary-buffer-koden är ett ofullständigt utkast. Slutför och granska den innan körning. Den första ändrade kontrollen är den redan motiverade förlängningen av numeriska högerranden från `1.2` till `1.8`, med bevarade ursprungliga källkolumner: högst två referenskörningar och 60 CPU-sekunder. Därefter verifieras en icke-noll historia och dess tangent.

Om tidsfelet därefter dominerar används den publicerade fjärde ordningens commutator-free Magnus-metoden på samma utökade system, med egen kontroll av fel och stabilitet på denna operator. Metodens formella ordning räknas inte som ett felcertifikat. [Blanes–Moan, ekvation 43](https://personales.upv.es/~serblaza/2006APNUM.pdf).

Bindningar ska täcka faktisk provtagen geometri, profil, initialfält, källidentitet och kanalparametrar. En ändrad profil får inte återanvända en gammal state-cache bara för att amplituderna är oförändrade.

**Leverans:** ny tillståndsägare, fokuserade tester, numeriskt register, sparade kontrollarrayer och verifierare. Ingen drift subtraheras som fysisk kompensation.

### Fas B — Slut den felbudget som den lokala gaten behöver

Bygg en gemensam källa för båda råa constrainterna:

\[
\mathcal E[g]
=\mathcal E_{\rm baseline}
+\Delta\mathcal E_{\rm local+reference}[g]
+\Delta G[C_\Sigma[g]].
\]

Baseline och geometribidrag läggs till en gång. Matter-skillnader summeras över den deklarerade inventeringen med ursprungliga vikter, tecken, degeneraciteter och koherenser.

Arbetet delas efter den verkliga återstående osäkerheten:

- Återanvänd certifierade bakgrundsregioner och komplettera eller begränsa återstående low/subgap-bidrag.
- Beräkna eller begränsa den ändrade historiens bidrag för samtliga berörda familjer. De tre \(\ell=0\)-gruppernas exakt noll state-respons i denna radiusfamilj återanvänds analytiskt; deras baseline finns kvar.
- Koppla det evolverade tillståndets högenergibeteende till den ägda subtraktionen och härled dess återstående svansgräns. Bakgrundens svanscertifikat överförs inte automatiskt till den nya historien.
- För vidare prepareringsfel, rums- och tidsfel, axialderivatans fel, kvadratur, koefficientintervall och avrundning till slutresidualen på \(I\).

Börja med gränser som kan avgöra gaten. Högre precision räknas fram endast när den kan ändra beslutet. Felbudgeten ska täcka den deklarerade källan och det lokala påståendet; den ska inte växa till ett krav på global matching eller full observationsanpassning.

**Leverans:** ett komplett register över täckning och fel för den nya tillståndslagen på \(I\). Saknade gränser förblir uttryckliga och blockerar certifiering.

### Fas C — Slut den bekräftade lokala incoming-gaten

Sök över två sammanhängande funktioner \(w(z)\) och \(U(z)\), genom den befintliga familjen

\[
\delta r=\chi(s)\left[s\,w(z)+\frac{s^3}{6}U(z)\right],
\qquad s=T-T_\Sigma.
\]

N, beta, a, källpreparation, intrinsisk incoming-identifikation och den fysiska seamens läge behålls. Historien använder den ägda retarderade prepareringsutvidgningen och dess kausala buffert. Constraintpåståendet gäller \(I\).

Som första numeriska domän används den befintliga kontrollens koordinater: ett inre intervall \(I=S(1)+[0.12,0.18]\), med slät axial övergång till noll inom \(S(1)+[0.09,0.21]\). Dessa är deklarerade koordinat- och beräkningsdomäner, ingen kosmologisk duration.

Representera funktionerna med Chebyshev-baser gånger en gemensam slät cutoff. Samtliga normal- och blandderivator kommer från dessa funktioner. Börja med åtta koefficienter per funktion; gå till 16 och 32 endast när residualens upplösningsfel motiverar det.

Lös båda ekvationerna samtidigt med en dämpad Newton-/Gauss–Newton-metod:

- Den fulla Jacobianen innehåller den retarderade state-responsen.
- Den certifierade geometriska principaloperatorn och dess fluxform används som preconditioner.
- Varje accepterat steg bevarar positiv radius och den deklarerade prepareringsdomänen.
- Startgissningar är solverkontroller. Det slutliga \(g\) verifieras genom ny evolution från samma \(C_{\rm up}\).

Kontrollera residualen mellan lösningsnoderna med en begränsad interpolations-/derivatarest. En nodvis träff räcker inte.

Om sökningen stannar analyseras den faktiska orsaken: numeriskt fel, otillräcklig representation eller en nödvändig matematisk begränsning. Ett NON-EXISTENCE-certifikat måste avse den fördeklarerade lokala klassen på \(I\) och använda den evolverade staten. Att en opåverkad region utanför \(I\) fortfarande har sin gamla residual avslutar inte denna gate.

**Leverans:** antingen lokal EXISTENCE med ett verifierat \((g,C_\Sigma[g])\) och båda residualerna inom tolerans över \(I\), eller scoped lokal NON-EXISTENCE med ett nytt nödvändigt hinder. OPEN är ett dokumenterat mellanresultat och uppfyller inte Done.

### Fas D — Skriv den citerbara, sakliga artikeln

Utgå från det befintliga engelska manuskriptet och dess PDF-byggare. Ge artikeln följande läsordning:

1. Den konkreta frågan, den deklarerade lokala modellen och gate-resultatet.
2. Importerade GR-, horizon-, Dirac- och CTP-byggdelar med tillämpningsvillkor.
3. Nested-gradient/separated-zone-tolkningen.
4. Den omkopplade tillståndslagen, beräkningsmetoden och felbudgeten.
5. Det lokala resultatet, reproduktion och exakt räckvidd.

Återanvänd attributionen till exempelvis [Klichs determinantidentitet](https://arxiv.org/abs/cond-mat/0209642), [Gérard–Häfner–Wrochnas tillståndskonstruktion](https://arxiv.org/abs/2008.10995) och [Bronnikov–Dehnen–Melnikovs black-universe-geometrier](https://arxiv.org/abs/gr-qc/0611022). Deras ursprungliga domäner och källor ska framgå; tillämpningen på NSC kräver den koppling som faktiskt visas.

Redovisa vad som är fixerat, löst respektive antaget, inklusive fria historiefunktioner och randdata. Påståendet om färre fria tillsatser ska kunna granskas i denna redovisning.

ΛCDM används som observationsspråk och möjlig framtida jämförelse. Gate-resultatet bedöms genom den deklarerade actionen. Metaforer om antimateria, entropi och intighet får inte bära resultatargumentet utan motsvarande härledd relation. Koordinatgränser och invariant singularitet skiljs åt där texten berör dem.

**Abstract, README och release ska alla ange samma lokala resultat, intervall \(I\), historikklass och tillståndslag.** Global parent/child-matching, extended stationarity och metric timestep beskrivs som efterföljande arbete, utan att göra det lokala resultatets publicering beroende av dem.

**Leverans:** uppdaterat manus, citerings- och claimkopplingar, resultattabell, relevanta figurer från sparade data samt byggd PDF.

### Fas E — Publicera hela beviskedjan för det lokala resultatet

Använd den befintliga offentliga checkouten. Importera den slutliga labbcommitens gate-register och hela deras beroendekedja från Git-objekt, med hashkontroll.

Utöka publiceringskontrollerna för de nya registerformaten och payloadbeskrivningarna. Vetenskapliga originalregister ändras inte för att passa importverktyget. Det frysta historiska manifestet behålls som historiskt.

Använd `0.27.0` som nästa föreslagna publikationsversion och synkronisera aktuella versionsfält. Historiska releasetester ska kontrollera sina historiska underlag; nya tester kontrollerar den aktuella publikationen.

Kör publiceringskontroller, fokuserade tester, gate-verifieraren och deterministisk PDF-kontroll. Byggkontrollen ska verifiera både upprepade PDF-byten och deras bindning till manus och vetenskapliga indata. CI ska återanvända låsta register och köra den nya relevanta grinden; en dokumentationsändring ska inte starta om gamla vetenskapliga kampanjer.

Publicera därefter normalt till GitHub `main`, med release och PDF. Verifiera publicerad commit och PDF-hash. Publiceringsposten anger både labb-SHA och offentlig SHA samt att resultatet är en **lokal incoming-gate**.

**Leverans:** offentlig, reproducerbar release. Zenodo ingår inte som standard eftersom det är valfritt. Ingen arXiv-publicering ingår.

## 3. Gränssnitt, verifiering och bindande slutkriterier

Produktionsvägen får tre tydliga resultatobjekt:

| Objekt | Obligatoriskt innehåll |
|---|---|
| Prepared incoming state | \(F,\delta F,F_z,\delta F_z\), käll- och historikbindning samt spektraltäckning. |
| Incoming constraints | Båda residualerna på \(I\), full historikderivata och uppdelad felbudget. |
| Local incoming gate certificate | Deklarerad lokal klass, sammanhängande intervall \(I\) med positiv längd, state-lag, reproducerbar historia eller nödvändigt hinder, tolerans, felgränser och provenance. |

Fryst `C0`, godtycklig matter-dictionary eller ofullständig felbudget får inte ge ett fysiskt PASS.

För **lokal numerisk EXISTENCE** krävs komponentvis

\[
\sup_{z\in I}
\left(
|\widehat{\mathcal E}_B(z)|+\epsilon_B(z)
\right)
\le 3\times10^{-11},
\qquad B=N,\beta.
\]

För **scoped lokal NON-EXISTENCE** krävs en nödvändig relation med en strikt, felkontrollerad separation från noll över hela den angivna lokala historikklassen. Misslyckad optimering är inte ett sådant resultat.

Den integrerade verifieringen omfattar:

- Oförändrad upstream-källa och korrekt retarderad variation.
- Finita differenskontroller av den nya sammanhängande utvärderingen.
- Kontroller som upptäcker borttagen \(\delta C\), tappad koherens, dubbelräknade vikter och felaktig substitution \(k=-E\).
- Numerisk rand-, fält- och derivatakontroll med tydligt skilda indikatorer och felgränser.
- Full deklarerad källtäckning och residualkontroll mellan noder på \(I\).
- Replay av slutcertifikatet samt publiceringens data-, manus- och PDF-bindningar.
- Samma uttryckliga lokala räckvidd i certifikat, abstract, README och release.

Varje verklig leverans verifieras och committas lokalt med residualer, toleranser, PASS/FAIL/OPEN och nästa namngivna lucka. Grok används för avgränsade delar och oberoende granskning; integrationen och slutpåståendet hålls samman.

**Sprinten är Done och ska stoppas när följande finns:**

1. Omkopplad incoming-state från samma upstream-källa, med matter och \(\delta C_\Sigma[g]\) i båda constrainterna.
2. Ett nytt stängt lokalt gate-resultat under denna lag.
3. Citerbart manus och verifierad PDF som uttryckligen presenterar resultatet som lokalt.
4. Verifierad offentlig GitHub-version med hela resultatets beroendekedja.

A/q/Ω/ζ/V_full förblir låsta, seeds förblir kontroller och ingen fabricerad stress eller ny `Γ_rest` införs. **Global matching, extended stationarity, metric timestep, full ΛCDM-anpassning och prövning av angränsande teorier får inte läggas till som villkor för Done.**

## 4. Troliga hinder och planerad hantering

| Hinder | Hantering |
|---|---|
| Rand- och tidsfel är större än den relevanta state-responsen. | Slut den identifierade numeriska luckan före lösningssökningen. Behåll källan och mät det råa felet. |
| Low/subgap eller den ändrade historiens UV-svans dominerar felbudgeten. | Använd riktade analytiska gränser och validerade paneler för just dessa bidrag. Utöka inte precisionen över redan certifierade regioner. |
| State-beroendet gör problemet illakonditionerat eller icke-lokalt. | Använd full retarderad Jacobian, geometrisk preconditioner, fluxform och dämpade steg. Geometrisk principalinvertibilitet räknas inte ensam som existens i den kopplade klassen. |
| Ett lokalt resultat glider över till ett globalt påstående eller får nya globala slutkrav. | Behåll den bekräftade lokala avgränsningen i certifikat, abstract, README och release. Globalt arbete ligger efter sprinten. |
| Publiceringsformat och versionshistorik blockerar ett annars färdigt resultat. | Anpassa den befintliga snapshot- och verifieringskedjan och bevara historiska bytes. Bygg en sammanhängande aktuell release. |

Det vetenskapliga utfallet är ännu inte känt. Avslutet är däremot fastställt: **omkopplad incoming-state + ett stängt ärligt lokalt gate-resultat + citerbar PDF + verifierad offentlig beviskedja.**
