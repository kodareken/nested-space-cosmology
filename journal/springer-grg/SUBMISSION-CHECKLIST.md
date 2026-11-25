# Step-by-Step Submission Instructions (Springer *GRG* -> arXiv Unlock)

Everything is built, verified, and ready in `journal/springer-grg/`. Nothing has been sent or submitted externally without your action. Follow these exact steps when you are ready:

---

### Step 1 — Fill Your City (and Optional ORCID) in `sn-article.tex` (1 minute)
1. Open `journal/springer-grg/sn-article.tex` at line 36:
   ```latex
   \affil*[1]{\orgname{Independent Researcher}, \orgaddress{\city{YOUR_CITY}, \country{Sweden}}}
   ```
   *(If you don't have a free ORCID yet, register in 60 seconds at `https://orcid.org/register` and add `\orcid{0000-...}` right after `\sur{Ek}` on line 34 — Springer journals strongly like seeing a verified ORCID for independent authors!)*
2. Recompile `sn-article.pdf` if you edited line 34–36 (or ask me here and I will recompile it in 3 seconds).

---

### Step 2 — Send the Short Courtesy Update to Prof. Paul Tod (Optional, Recommended)
- Open `journal/springer-grg/NOTE-TO-PAUL-TOD.md`, attach `journal/springer-grg/sn-article.pdf`, and send to `paul.tod@sjc.ox.ac.uk`.
- Why do this right before or right as you submit to *GRG*? You promised him an update, *GRG* is one of his home journals, and keeping goodwill with an Oxford relativist who already endorsed you is priceless.

---

### Step 3 — Upload to Springer *General Relativity and Gravitation* (5 minutes, $0 Fee)
1. Go to the official *General Relativity and Gravitation* submission portal:  
   `https://link.springer.com/journal/10714` -> click **"Submit manuscript"** (Springer Nature Snapp).
2. Choose **Article Type:** `Research` (Original Research Article).
3. Upload files from `journal/springer-grg/`:
   - **Main Manuscript File:** `sn-article.tex`
   - **Support / Style / Bib Files:** `sn-jnl.cls`, `sn-mathphys-num.bst`, `sn-bibliography.bib` (and `sn-article.bbl` if requested)
   - **Figure Files:** `figures/fig1-bounce-and-geometry.pdf`, `figures/fig2-variational-budget-surface.pdf`, `figures/fig3-neck-cubic-and-cascade.pdf` *(600-DPI `.png` versions are also in `figures/` if the portal asks for raster previews)*
   - **Cover Letter:** Paste or upload `COVER-LETTER-GRG.md` (with your city filled in).
4. **Publishing Model Question (Important!):** When asked between Open Access (APC fee) vs. **Standard Subscription / Hybrid ($0 fee)**, select **Subscription ($0)** unless you specifically wish to pay for Open Access. Subscription articles in *GRG* receive a full Springer DOI (`10.1007/s10714-...`) and Springer's self-archiving policy explicitly allows posting the accepted author manuscript on **arXiv** immediately!

---

### Step 4 — Unlocking arXiv (`MOD-115301`) After Journal Acceptance
- Once *GRG* (or backup #2: *Classical and Quantum Gravity* / *European Physical Journal C*, which uses the **exact same `sn-jnl` files**) accepts the paper and assigns your `10.1007/...` DOI:
  1. Submit `sn-article.tex` + figures to `gr-qc` on arXiv, entering the published **Journal Reference** and **DOI** in the arXiv metadata fields.
  2. That satisfies the `MOD-115301` condition cleanly and establishes your peer-reviewed track record on `gr-qc` for future papers (including Paper II!).
