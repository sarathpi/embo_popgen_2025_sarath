library(admixtools)
library(tidypopgen)
f2_blocks = f2_from_precomp("./data/f2_tidypopgen", verbose = FALSE)

neand_euras_wave <- qpwave(data = f2_blocks,
      left = c("French","Spanish","Tuscan"),
      right = c("AltaiNea","Mota", "Yoruba", "Denisova", "Mbuti")
)
neand_euras_wave

french_adm <- qpadm(data = f2_blocks,
      left = c("Loschbour", "LBK", "Yamnaya"),
      right = c("Mbuti", "Mota", "Dinka", "Yoruba", "Han"),
      target= "French")

french_adm$popdrop

base_edges <- matrix(
  c("R",	"Mbuti",
    "R", "eAfr",
    "eAfr",	"Dinka",
    "eAfr",	"outAfrica",
    "outAfrica",	"Han",
    "outAfrica",	"Loschbour"),
  ncol=2,byrow = TRUE,
  dimnames=list(NULL, c("from","to")))

base_edges %>% edges_to_igraph() %>% plot_graph()

base_edges <- matrix(
  c("R",	"Mbuti",
    "R", "eAfr",
    "eAfr",	"Dinka",
    "eAfr",	"outAfrica",
    "outAfrica",	"Han",
    "outAfrica",	"Loschbour"),
  ncol=2,
  byrow = TRUE,
  dimnames=list(NULL, c("from","to")))

base_edges

base_igraph <- base_edges %>% edges_to_igraph()

is_valid(base_igraph)

base_igraph %>% plot_graph()

base_igraph %>% plotly_graph()

base_qpgraph <- qpgraph(data = f2_blocks, graph = base_igraph)

base_qpgraph$f3

base_qpgraph$f3 %>% filter(abs(z)>2)

base_qpgraph$edges %>% plot_graph()

fits = qpgraph_resample_multi(f2_blocks, 
                              graphlist = list(base_qpgraph[[1]], base_swapped_qpgraph[[1]]), 
                              nboot = 100)
compare_fits(fits[[1]]$score_test, fits[[2]]$score_test)

base_igraph %>% plot_graph(highlight_unidentifiable = TRUE)

yamnaya_edges <- matrix(
  c(
    "R",	"Mbuti",
    "R",	"eAfr",
    "eAfr", "Dinka",
    "eAfr", "outAfrica",
    "outAfrica", "Han",
    "outAfrica", "wEurasian",
    "wEurasian", "Yamnaya",
    "wEurasian", "Loschbour"),
  ncol = 2,
  byrow = TRUE,
  dimnames = list(NULL, c("from", "to")))
yamnaya_igraph <- yamnaya_edges %>% edges_to_igraph()
yamnaya_igraph %>% plot_graph()

lbk_extra_edges <- matrix(
  c(
    "R",	"Mbuti",
    "R",	"eAfr",
    "eAfr", "pBasalEurasian",
    "eAfr", "Dinka",
    "pBasalEurasian", "BasalEurasian",
    "pBasalEurasian","outAfrica",
    "outAfrica", "Han",
    "outAfrica","wEurasian",
    "wEurasian", "Yamnaya",
    "wEurasian", "pLoschbour",
    "pLoschbour", "Loschbour",
    "pLoschbour","WHG",
    "BasalEurasian", "pLBK",
    "WHG", "pLBK",
    "pLBK","LBK"),
  ncol = 2,
  byrow = TRUE,
  dimnames = list(NULL, c("from", "to")))
lbk_extra_igraph <- lbk_extra_edges %>% edges_to_igraph()
lbk_extra_igraph %>% plot_graph()

is_valid(lbk_extra_igraph)

lbk_extra_igraph %>% plot_graph(highlight_unidentifiable = TRUE)

lbk_extra_qpgraph <- qpgraph(data = f2_blocks, graph = lbk_extra_igraph)
lbk_extra_qpgraph$edges %>% plot_graph()

