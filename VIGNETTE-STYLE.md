# Vignette conventions

This file records the conventions the dkge vignette suite already follows, and
the rules agreed for the ones it does not. Every number below was measured at
commit 276adcc and reproduced by two people.

Read this before editing a vignette. It exists so that the next person does not
rediscover a convention as a defect, which is what happened three times while
these rules were being written.

## Page structure

Every article fills four slots in order.

1. What this page is for, and what you must already have read.
2. The worked example.
3. What you may and may not claim from it.
4. Where to go next.

Slot 1 is one paragraph. It opens on the reader's situation, not on a
description of the document. Prerequisites follow in the same paragraph. The
scope contract belongs in the `description:` field of the YAML header, where
pkgdown renders it in the article index.

Slot 4 is named "Where to go next" and it replaces the closing section. Do not
append it below a `## Summary`. A recap section has no grammatical slot for
pointing anywhere, which is why every recap-closing article in the suite ends
without an onward link while linking freely in its body.

## Links

Use `vignette("dkge-partial-effect-spaces")` everywhere, including the
"Where to go next" section. Do not write raw `.html` links.

Both measurements behind that rule, because a rule without its evidence gets
reversed by the next person who thinks about it, and this one was reversed
twice:

- **On the pkgdown site**, downlit autolinks the call. A one-article build of
  `dkge-workflow` produced five anchors of the form
  `<a href="../articles/dkge.html">vignette("dkge")</a>` and zero inert code
  spans. It is a real link.
- **In the installed doc directory**, there is no downlit, so it renders as
  code. A reader who types it gets the vignette, so the reference still works —
  as a command rather than a hyperlink.

A relative `.html` link is clickable in both places, which is why this looked
like the better choice for a while. It loses on three counts: it is dead when
someone reads the `.Rmd` source on GitHub, where nothing is built; it is not
what an R user types; and dplyr, tibble, knitr and rmarkdown carry zero raw
`.html` links in vignette source.

Every article's final section carries at least one `vignette()` reference, and
no raw `.html` links.

## Voice

Address the reader. Fourteen of the eighteen articles already do.

Use `we` only for what the authors built for the example: "We simulate 10
subjects". Keep it to a handful of instances and keep it out of instruction.

Avoid the lecturer register: "we will proceed through four main steps", "let us
consider", "this will allow us to demonstrate". Thirty-five such markers survive
in four files and none appear in the other fourteen.

Put the identifier in subject position. "First, DKGE constructs a projector
through `dkge_projector_K()`" buries the actor, so the sentence needs "First" to
place itself. Start with `dkge_projector_K()` and the ordering explains itself.
Sequence words fall away once the subject is right.

## Explaining code

Explain before you show. Sixteen of the eighteen articles put prose above the
chunk that produces output, and that is the house convention.

Reader guidance belongs in prose, not in code comments. A comment that tells the
reader what to do in a real analysis is invisible to anyone skimming and to
every screen reader.

Every figure carries `fig.alt`. All 32 currently do. A figure also needs prose
that a reader arriving at the image can use.

## Terminology

Define a term where a reader first meets it. The debut is the first occurrence
outside backticks and outside headings; an identifier is not a use of the word.
A definition may live in a heading.

`salience`, `estimand`, and `q-space` are used across the suite and defined
nowhere. `design kernel` is defined three articles after its debut.

## Naming

One name per role, across all eighteen articles: `toy` for simulated input,
`bundle` for the object `dkge_data()` returns, `fit` for the fitted model.

Where a page holds two fits, qualify both names. A bare `fit` beside a
`fit_analytic` invites prose to drift onto the wrong object, which is how two
sentences in `dkge-partial-effect-spaces` came to describe an experiment the
reader was not looking at.

## Hedging

State the positive claim narrowly enough that a denial is unnecessary.

Keep an explicit denial where a competent reader would otherwise reach the wrong
reading. Simulation-recovery figures qualify. Most of the fifty-six instances in
the suite do not, and a denial the reader was never going to need still plants
the reading it rejects.

## Spelling

The package uses en-US in prose. Match argument names to the API they document:
an argument called `normalize` is described as normalizing.

Some British spellings are ggplot2 argument names and stay as they are.

## Checks

`inst/checks/` holds one extractor and several rules. The extractor is the
contract: it strips the YAML header, fenced chunks, table rows, raw HTML, inline
code, quoted spans, and blockquotes before any rule matches. Without the last
three, the first document a rule flags is the one that defines it.

Three rules gate: bare figures, register markers, opener shape. Three report a
number without gating: closing-section type, line width, denial density. One
warns only: reader prose inside code comments.

Match terms by lemma. A rule written for `cross-fitting` misses the article that
says `cross-fitted`.

## Changing these conventions

Additive rules travel on their own. The `description:` field reached eighteen
of nineteen articles with no enforcement, because it adds something unambiguous
and costs one line.

Substitutive rules do not. A rule that asks people to stop doing something, or
to prefer one form over another that looks equally reasonable, needs a check
behind it. Two such rules were reversed mid-work here and both times someone
kept executing the old version for an hour, because recording a reversal is not
the same as delivering it.

So: if a convention in this file replaces an existing practice, it ships with a
check. If it only adds, a sentence here is enough.
