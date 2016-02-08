require(git2r)
require(magrittr)
repo <- repository(".")
sha <- head(repo) %>% branch_target
cat(file = "VERSION.tex", sep="",
  "\\providecommand{\\VERSION}{",
  lookup(repo, sha) %>% when %>% strtrim(10),
  ", commit ",
  strtrim(sha, 7),
  "}\n")
