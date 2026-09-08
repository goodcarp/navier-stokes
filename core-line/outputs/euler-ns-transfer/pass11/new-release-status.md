# Public release status, 8 September 2026

Checked on **2026-09-08 at approximately 13:20–13:23 UTC (09:20–09:23 EDT)**. This is a bounded check of primary release channels, not an audit of the Euler proof and not evidence about unpublished work.

**No newly public ordinary or hypodissipative Navier–Stokes manuscript was located.** The current author-hosted Euler PDF and repository main commit are unchanged from the artifacts already inspected in this project.

| Primary artifact/channel | Fresh observation |
|---|---|
| [Author-hosted Euler manuscript](https://cims.nyu.edu/~tristanb/euler.pdf) | Successfully downloaded directly. It is the same 112-page *Blowup for the Euler Equations with Smooth Forcing*, byte-for-byte identical to the archived copy. Its abstract concerns forced incompressible Euler on R3, with the force smooth through blowup; it does not state an ordinary-NS theorem. |
| [Public fluid_lean repository](https://github.com/tristanbuckmaster/fluid_lean) and [commit history](https://github.com/tristanbuckmaster/fluid_lean/commits/main/) | The current main commit, independently fetched through the public GitHub API, remains `d0124689230b58b4f86e7b90ac59de06404b3b6b`, dated 2026-09-08 04:07:52 UTC. It has no parent commit. The visible root contains `affinecore`, `boussinesq-blowup`, and `euler-blowup`; no new NS project is listed. |
| [Buckmaster publication page](https://cims.nyu.edu/~tristanb/publications/) and [Alpöge homepage](https://people.math.harvard.edu/~alpoge/) | No newly posted ordinary/hypodissipative NS paper was found on the returned pages. These publication lists are not themselves reliable evidence of the absence of a new release: the already available Euler PDF is not surfaced in the publication list either. |
| [arXiv author search: Alpoge](https://arxiv.org/search/?query=Alpoge&searchtype=author&abstracts=show&order=-announced_date_first&size=50) | Directly fetched 21 returned records, including namesakes. No Euler, Boussinesq, porous-media, or NS title appeared. |
| [arXiv author search: Buckmaster](https://arxiv.org/search/?query=Buckmaster&searchtype=author&abstracts=show&order=-announced_date_first&size=50) | Directly fetched 34 returned records, including namesakes. The NS titles returned were the older weak-solution/nonuniqueness papers; no new smooth-data ordinary or hypodissipative NS manuscript appeared. |

Targeted searches restricted to arXiv, the author domains, the release repository, and official OpenAI release pages likewise located no new ordinary-NS proof artifact. Community posts and secondary reports were not used as evidence that a proof exists. Search-index and release-timing gaps remain possible.

The hypodissipative statement in the user's supplied announcement is a report of ongoing work, explicitly described there as unreleased with unfinished Lean verification. It should retain that status unless an actual manuscript or proof artifact becomes publicly inspectable. The stronger ordinary smooth-forced NS report in the second screenshot is likewise not a public theorem established by the artifacts found in this check. No conclusion about the authors' unpublished progress follows.

For competitive planning, this check supplies no reason to treat an ordinary smooth-forced NS target as already publicly completed by this release. It also supplies no basis for predicting whether our project can finish first. Smooth forcing remains an allowed target when the precise theorem's data, regularity, domain, and energy requirements are met; the current unforced experiment does not impose a ban on that route.

Provenance retained in the original collection workspace (the packaged pass
publishes URLs, hashes and selected JSON records without duplicating the
downloaded PDF and HTML):

- `euler-current.pdf`: 1,105,306 bytes; SHA256 `97ef408bff09b4f6ed9f3867734d1eb2245f3f34e6334b28136c84c02d0ae8d8`, identical to `work/sources/euler.pdf`.
- `fluid-current-commit.json`: fresh response from the public main-commit API; tree SHA `6a5e38b2932120637b2f854251d1a09d3eedeaad`.
- `arxiv-alpoge-current.html` and `arxiv-buckmaster-current.html`, with compact parsed title/date/link inventories in the corresponding `*-results.json` files.

Direct web extraction of the PDF and arXiv search pages failed, so the successful fresh public downloads were used for the file comparison and author-list checks. No private source, message, or account was accessed, and no contact or recurring monitoring was initiated.
