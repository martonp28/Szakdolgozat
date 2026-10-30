library(dplyr)
library(tidyr)
library(network)
library(ergm)
library(ergm.multi)
library(tictoc)
library(writexl)
library(readxl)
library(sna)

setwd("C:/ELTE/Szakdoga/data/szakdoga ml ergm")

###########################################################################
# Beolvasása az adatoknak

# Trade
trade_sheets <- excel_sheets("networks.xlsx")
trade <- list()

for (sheet in trade_sheets) {
  
  df <- read_excel("networks.xlsx", sheet = sheet)
  
  stopifnot("Sender" %in% names(df))
  
  sender <- as.character(df$Sender)
  country_cols <- setdiff(names(df), "Sender")
  
  # Ellenőrizzük, hogy az oszlopok tényleges országnevei
  # megegyeznek-e a kibocsátó országokkal.
  stopifnot(
    !anyDuplicated(sender),
    !anyDuplicated(country_cols),
    setequal(country_cols, sender)
  )
  
  mx <- as.matrix(df[, country_cols, drop = FALSE])
  rownames(mx) <- sender
  
  # Az oszlopokat a sorok sorrendjére rendezzük,
  # nem egyszerűen átnevezzük.
  mx <- mx[, sender, drop = FALSE]
  storage.mode(mx) <- "numeric"
  
  stopifnot(
    identical(rownames(mx), colnames(mx)),
    all(is.finite(mx)),
    all(mx >= 0)
  )
  
  trade[[sheet]] <- mx
}
names(trade)
trade[['2015']]
dim(trade[["2015"]])
head(trade[["2015"]])
rownames(trade[["2015"]])

# World Bank
wb_sheets <- excel_sheets('worldbank_annual.xlsx')
wb <- list()

for (sheet in wb_sheets){
  df <- read_excel(
    "worldbank_annual.xlsx",
    sheet = sheet
  )
  
  df <- as.data.frame(df)
  
  rownames(df) <- df$Node
  df$Node <- NULL
  
  wb[[sheet]] <- df
}
wb[["2016"]]
names(wb)
dim(wb[["2015"]])
rownames(wb[["2015"]])

# CEPII
cepii_sheets <- excel_sheets("cepii.xlsx")
cepii <- list()

for (sheet in cepii_sheets) {
  
  df <- read_excel("cepii.xlsx", sheet = sheet)
  
  stopifnot("Sender" %in% names(df))
  
  sender <- as.character(df$Sender)
  country_cols <- setdiff(names(df), "Sender")
  
  stopifnot(
    !anyDuplicated(sender),
    !anyDuplicated(country_cols),
    setequal(country_cols, sender)
  )
  
  mx <- as.matrix(df[, country_cols, drop = FALSE])
  rownames(mx) <- sender
  
  mx <- mx[, sender, drop = FALSE]
  storage.mode(mx) <- "numeric"
  
  stopifnot(
    identical(rownames(mx), colnames(mx)),
    all(is.finite(mx))
  )
  
  cepii[[sheet]] <- mx
}
cepii[["log_dist"]]
names(cepii)
dim(cepii[["log_dist"]])
rownames(cepii[["log_dist"]])


###############################################################################
# MULTILAYER ERGM - 2015-2024
# Teljes Model 8 specifikáció éves ciklusban
#
# Feltétel:
# trade, wb, cepii, routes már be vannak töltve.
###############################################################################

# setwd("C:/ELTE/Szakdoga/data/ml ergm/nested")

###############################################################################
# 0. BEÁLLÍTÁSOK
###############################################################################

names(wb[["2016"]])

dim(trade[["2015"]])
dim(wb[["2015"]])

dim(cepii[["log_dist"]])
dim(cepii[["comlang"]])
dim(cepii[["colony"]])
dim(cepii[["contig"]])

# sapply(routes, dim)

identical(
  rownames(trade[["2015"]]),
  rownames(wb[["2015"]])
)

identical(
  rownames(trade[["2015"]]),
  rownames(cepii[["log_dist"]])
)

years <- 2015:2024
n_nodes <- 36
n_layers <- 4
n_layer_gof_sim <- 500L

# years <- 2024L
# n_layer_gof_sim <- 50L

# Kimeneti mappák
dir.create("model diagnostics", showWarnings = FALSE)
dir.create("annual_models", showWarnings = FALSE)
dir.create("annual_gof", showWarnings = FALSE)
dir.create("annual_layers", showWarnings = FALSE)
dir.create("annual_logs", showWarnings = FALSE)
dir.create("annual_results", showWarnings = FALSE)

###############################################################################
# 1. STATIKUS NODE-SORREND ÉS CEPII / ROUTE VÁLTOZÓK
###############################################################################

# Node-sorrend a trade hálózatból
node_order <- rownames(trade[[as.character(years[1])]])

# megnézi hogy pontosan annyi ország van-e a node orderben amennyit várunk
stopifnot(
  length(node_order) == n_nodes,
  identical(node_order, colnames(trade[[as.character(years[1])]]))
)
# nincs output -> jó

# Minden évben ugyanaz legyen a trade node-sorrend
for (year in years) {
  y <- as.character(year)
  
  stopifnot(
    identical(rownames(trade[[y]]), node_order),
    identical(colnames(trade[[y]]), node_order)
  )
}

# WB változók ellenőrzése
for (year in years) {
  
  wb_year <- wb[[as.character(year)]]
  
  stopifnot(
    all(c("log_hightech", "log_gdp", "log_labor") %in% names(wb_year)),
    identical(rownames(wb_year), node_order)
  )
  
  stopifnot(
    !anyNA(wb_year$log_hightech),
    !anyNA(wb_year$log_gdp),
    !anyNA(wb_year$log_labor)
  )
}

###############################################################################
# 2. STATIKUS CEPII VÁLTOZÓK
###############################################################################

comlang <- cepii[["comlang"]]
log_dist <- cepii[["log_dist"]]
colony <- cepii[["colony"]]
contig <- cepii[["contig"]]

# CEPII node-sorrend
stopifnot(
  identical(rownames(comlang), colnames(comlang)),
  identical(rownames(log_dist), colnames(log_dist)),
  identical(rownames(colony), colnames(colony)),
  identical(rownames(contig), colnames(contig)),
  
  identical(rownames(comlang), node_order),
  identical(rownames(log_dist), node_order),
  identical(rownames(colony), node_order),
  identical(rownames(contig), node_order)
)

# Hiányzó értékek
stopifnot(
  !anyNA(comlang),
  !anyNA(log_dist),
  !anyNA(colony),
  !anyNA(contig)
)

###############################################################################
# 4. SEGÉDFÜGGVÉNYEK
###############################################################################

# 36x36 -> 144x144 blokkdiagonális mátrix
make_block144 <- function(mx, nodes, n_layers = 4) {
  
  n <- length(nodes)
  
  stopifnot(
    nrow(mx) == n,
    ncol(mx) == n
  )
  
  out <- matrix(
    0,
    nrow = n * n_layers,
    ncol = n * n_layers
  )
  
  for (l in 0:(n_layers - 1)) {
    
    idx <- (1:n) + l * n
    
    out[idx, idx] <- mx
  }
  
  layer_nodes <- unlist(
    lapply(
      1:n_layers,
      function(l) paste0("L", l, "_", nodes)
    )
  )
  
  dimnames(out) <- list(layer_nodes, layer_nodes)
  
  out
}


# Mutual dyadok száma egy 36x36 directed rétegben
mutual_dyads <- function(mx) {
  sum(mx == 1 & t(mx) == 1) / 2
}

# A multilayer szimulációból négy 36x36-os réteg kinyerése
extract_layer_matrices <- function(
    net,
    n = 36L,
    n_layers = 4L
) {
  
  A <- as.matrix(
    net,
    matrix.type = "adjacency"
  )
  
  expected_dim <- c(n * n_layers, n * n_layers)
  
  if (!identical(dim(A), expected_dim)) {
    stop(
      "A szimulált hálózat dimenziója nem 144x144. ",
      "Ellenőrizd a Layer objektum reprezentációját."
    )
  }
  
  layers <- lapply(seq_len(n_layers), function(l) {
    
    idx <- ((l - 1L) * n + 1L):(l * n)
    
    1L * (A[idx, idx, drop = FALSE] != 0)
  })
  
  names(layers) <- paste0("L", seq_len(n_layers))
  
  layers
}


# Fokszámgyakoriságok és irányított triádcenzus egyetlen rétegben
layer_gof_stats <- function(A) {
  
  A <- 1L * (A != 0)
  
  n <- nrow(A)
  
  stopifnot(
    ncol(A) == n,
    all(diag(A) == 0)
  )
  
  in_freq <- tabulate(
    as.integer(colSums(A)) + 1L,
    nbins = n
  )
  
  out_freq <- tabulate(
    as.integer(rowSums(A)) + 1L,
    nbins = n
  )
  
  names(in_freq) <- paste0(
    "indegree_", 0:(n - 1L)
  )
  
  names(out_freq) <- paste0(
    "outdegree_", 0:(n - 1L)
  )
  
  net <- network::network(
    A,
    directed = TRUE,
    matrix.type = "adjacency",
    loops = FALSE
  )
  
  triads <- sna::triad.census(
    net,
    mode = "digraph"
  )
  
  triad_values <- if (is.matrix(triads)) {
    as.numeric(triads[1, ])
  } else {
    as.numeric(triads)
  }
  
  triad_names <- c(
    "003", "012", "102", "021D",
    "021U", "021C", "111D", "111U",
    "030T", "030C", "201", "120D",
    "120U", "120C", "210", "300"
  )
  
  if (length(triad_values) != 16L) {
    stop("A triadcensus nem 16 statisztikát adott vissza.")
  }
  
  names(triad_values) <- paste0(
    "triad_", triad_names
  )
  
  c(
    in_freq,
    out_freq,
    triad_values
  )
}


# Megfigyelt és szimulált rétegek összehasonlítása
compute_layer_gof <- function(
    observed_layers,
    sim_networks,
    year
) {
  
  if (inherits(sim_networks, "network")) {
    sim_networks <- list(sim_networks)
  }
  
  if (!is.list(sim_networks) || length(sim_networks) < 20L) {
    stop("Nem megfelelő számú szimulált hálózat érkezett a layer-GOF-ba.")
  }
  
  observed_stats <- lapply(
    observed_layers,
    layer_gof_stats
  )
  
  # Csak a statisztikákat tároljuk, nem a mátrixokat.
  simulated_stats <- lapply(sim_networks, function(net) {
    
    mats <- extract_layer_matrices(net)
    
    lapply(mats, layer_gof_stats)
  })
  
  result <- list()
  k <- 0L
  
  for (layer_name in names(observed_stats)) {
    
    obs <- observed_stats[[layer_name]]
    
    sim_matrix <- do.call(
      rbind,
      lapply(
        simulated_stats,
        function(x) x[[layer_name]]
      )
    )
    
    if (!identical(colnames(sim_matrix), names(obs))) {
      stop(
        "A megfigyelt és szimulált statisztikák nem egyeznek: ",
        layer_name
      )
    }
    
    for (stat_name in names(obs)) {
      
      observed_value <- as.numeric(obs[[stat_name]])
      sim_values <- as.numeric(sim_matrix[, stat_name])
      
      p_lower <- (
        sum(sim_values <= observed_value) + 1
      ) / (length(sim_values) + 1)
      
      p_upper <- (
        sum(sim_values >= observed_value) + 1
      ) / (length(sim_values) + 1)
      
      mc_p <- min(
        1,
        2 * min(p_lower, p_upper)
      )
      
      if (grepl("^indegree_", stat_name)) {
        
        statistic <- "indegree_frequency"
        category <- sub("^indegree_", "", stat_name)
        
      } else if (grepl("^outdegree_", stat_name)) {
        
        statistic <- "outdegree_frequency"
        category <- sub("^outdegree_", "", stat_name)
        
      } else {
        
        statistic <- "triad_census"
        category <- sub("^triad_", "", stat_name)
      }
      
      q <- quantile(
        sim_values,
        probs = c(0.025, 0.975),
        na.rm = TRUE,
        names = FALSE
      )
      
      k <- k + 1L
      
      result[[k]] <- data.frame(
        Year = year,
        Layer = layer_name,
        Statistic = statistic,
        Category = category,
        Observed = observed_value,
        Sim_Mean = mean(sim_values),
        Sim_Lower_95 = q[1],
        Sim_Upper_95 = q[2],
        MC_two_sided_p = mc_p,
        N_simulations = length(sim_values),
        stringsAsFactors = FALSE
      )
    }
  }
  
  dplyr::bind_rows(result)
}

# ERGM fit közben:
# - a konzolban látszik a progress
# - ugyanaz az output logfájlba kerül
# - a warningokat összegyűjtjük
# fit_ergm_logged <- function(formula, log_file) {
#   
#   warnings <- character()
#   err <- NULL
#   model <- NULL
#   
#   sink(log_file, split = TRUE)
#   
#   model <- tryCatch(
#     
#     withCallingHandlers(
#       
#       ergm(formula),
#       
#       warning = function(w) {
#         
#         warnings <<- c(
#           warnings,
#           conditionMessage(w)
#         )
#         
#         invokeRestart("muffleWarning")
#       }
#     ),
#     
#     error = function(e) {
#       err <<- e
#       NULL
#     },
#     
#     finally = {
#       sink()
#     }
#   )
#   
#   fit_output <- readLines(
#     log_file,
#     warn = FALSE
#   )
#   
#   if (!is.null(err)) {
#     stop(err)
#   }
#   
#   list(
#     model = model,
#     output = fit_output,
#     warnings = unique(warnings)
#   )
# }


# MCMLE outputból a végső convergence információk
parse_mcmle_info <- function(lines) {
  
  lines <- trimws(lines)
  
  # Konvergencia p-értékének keresése
  conv_idx <- grep(
    "Convergence test p-value\\s*:",
    lines,
    perl = TRUE
  )
  
  conv_p <- NA_real_
  conv_p_reported <- NA_character_
  
  if (length(conv_idx) > 0) {
    
    last_line <- lines[tail(conv_idx, 1)]
    
    p_raw <- sub(
      ".*Convergence test p-value\\s*:\\s*",
      "",
      last_line,
      perl = TRUE
    )
    
    p_match <- regmatches(
      p_raw,
      regexpr(
        "<\\s*[0-9]+(?:\\.[0-9]*)?(?:[eE][+-]?[0-9]+)?|[0-9]+(?:\\.[0-9]*)?(?:[eE][+-]?[0-9]+)?",
        p_raw,
        perl = TRUE
      )
    )
    
    if (length(p_match) > 0 && nzchar(p_match)) {
      
      conv_p_reported <- gsub(
        "\\s+",
        "",
        p_match
      )
      
      # Ha az output felső korlátot közöl, nem
      # állítjuk be pontos p-értékként.
      if (!grepl("^<", conv_p_reported)) {
        conv_p <- suppressWarnings(
          as.numeric(conv_p_reported)
        )
      }
    }
  }
  
  # A konvergenciastátusz külön sorban is szerepelhet.
  success <- any(
    grepl(
      "Converged with 99% confidence",
      lines,
      fixed = TRUE
    )
  )
  
  failure <- any(
    grepl(
      "Failed to converge|Did not converge",
      lines,
      ignore.case = TRUE
    )
  )
  
  if (success) {
    converged_99 <- TRUE
  } else if (failure) {
    converged_99 <- FALSE
  } else {
    converged_99 <- NA
  }
  
  # Iterációszám, a sor eleji szóközök megengedésével
  iteration_lines <- grep(
    "^\\s*Iteration\\s+[0-9]+",
    lines,
    value = TRUE,
    perl = TRUE
  )
  
  if (length(iteration_lines) == 0) {
    
    n_iter <- NA_integer_
    
  } else {
    
    iteration_numbers <- suppressWarnings(
      as.integer(
        sub(
          "^\\s*Iteration\\s+([0-9]+).*",
          "\\1",
          iteration_lines,
          perl = TRUE
        )
      )
    )
    
    iteration_numbers <- iteration_numbers[
      !is.na(iteration_numbers)
    ]
    
    n_iter <- if (length(iteration_numbers) > 0) {
      max(iteration_numbers)
    } else {
      NA_integer_
    }
  }
  
  list(
    MCMLE_convergence_p = conv_p,
    MCMLE_convergence_p_reported = conv_p_reported,
    MCMLE_converged_99 = converged_99,
    MCMLE_iterations = n_iter
  )
}


# summary(model) alapján fit mutatók
parse_model_fit <- function(model) {
  
  txt <- capture.output(
    summary(model)
  )
  
  # AIC/BIC/MC SE
  aic_line <- grep(
    "^\\s*AIC:",
    txt,
    value = TRUE
  )
  
  if (length(aic_line) == 0) {
    
    AIC_value <- NA_real_
    BIC_value <- NA_real_
    MC_SE <- NA_real_
    
  } else {
    
    AIC_value <- suppressWarnings(
      as.numeric(
        sub(
          ".*AIC:\\s*([0-9.eE+-]+).*",
          "\\1",
          aic_line[1]
        )
      )
    )
    
    BIC_value <- suppressWarnings(
      as.numeric(
        sub(
          ".*BIC:\\s*([0-9.eE+-]+).*",
          "\\1",
          aic_line[1]
        )
      )
    )
    
    MC_SE <- suppressWarnings(
      as.numeric(
        sub(
          ".*MC Std\\. Err\\. =\\s*([0-9.eE+-]+).*",
          "\\1",
          aic_line[1]
        )
      )
    )
  }
  
  # Null deviance
  null_line <- grep(
    "Null Deviance:",
    txt,
    value = TRUE
  )
  
  if (length(null_line) == 0) {
    
    Null_Deviance <- NA_real_
    Null_df <- NA_real_
    
  } else {
    
    Null_Deviance <- suppressWarnings(
      as.numeric(
        sub(
          ".*Null Deviance:\\s*([0-9.eE+-]+).*",
          "\\1",
          null_line[1]
        )
      )
    )
    
    Null_df <- suppressWarnings(
      as.numeric(
        sub(
          ".*on\\s+([0-9.eE+-]+)\\s+degrees.*",
          "\\1",
          null_line[1]
        )
      )
    )
  }
  
  # Residual deviance
  residual_line <- grep(
    "Residual Deviance:",
    txt,
    value = TRUE
  )
  
  if (length(residual_line) == 0) {
    
    Residual_Deviance <- NA_real_
    Residual_df <- NA_real_
    
  } else {
    
    Residual_Deviance <- suppressWarnings(
      as.numeric(
        sub(
          ".*Residual Deviance:\\s*([0-9.eE+-]+).*",
          "\\1",
          residual_line[1]
        )
      )
    )
    
    Residual_df <- suppressWarnings(
      as.numeric(
        sub(
          ".*on\\s+([0-9.eE+-]+)\\s+degrees.*",
          "\\1",
          residual_line[1]
        )
      )
    )
  }
  
  list(
    AIC = AIC_value,
    BIC = BIC_value,
    MC_SE = MC_SE,
    Null_Deviance = Null_Deviance,
    Null_df = Null_df,
    Residual_Deviance = Residual_Deviance,
    Residual_df = Residual_df
  )
}


# Coefficientek kivétele
extract_coefficients <- function(model, year) {
  
  cm <- coef(summary(model))
  
  out <- data.frame(
    Year = year,
    Term = rownames(cm),
    Estimate = cm[, 1],
    Std_Error = cm[, 2],
    MCMC_pct = cm[, 3],
    z_value = cm[, 4],
    p_value = cm[, 5],
    stringsAsFactors = FALSE
  )
  
  out
}


# Paraméterek csoportosítása
classify_term <- function(term) {
  
  if (grepl("^L\\(", term) && grepl("edges", term)) {
    return("Layer_edges")
  }
  
  if (grepl("^L\\(", term) && grepl("mutual", term)) {
    return("Layer_mutual")
  }
  
  if (grepl("^nodeo", term) || grepl("^nodei", term)) {
    return("Node_covariate")
  }
  
  if (grepl("^edgecov\\.", term)) {
    return("Dyadic_control")
  }
  
  return("Other")
}

# Szebb label a többi változónak
term_labels <- c(
  "nodeocov.log_hightech" = "High-tech export (outgoing)",
  "nodeicov.log_hightech" = "High-tech export (incoming)",
  "nodeocov.log_gdp" = "GDP per capita (outgoing)",
  "nodeicov.log_gdp" = "GDP per capita (incoming)"
)

# A fenti névduplikáció miatt külön kezeljük
term_label <- function(term) {
  
  if (term == "nodeocov.log_hightech") {
    return("High-tech export (outgoing)")
  }
  
  if (term == "nodeicov.log_hightech") {
    return("High-tech export (incoming)")
  }
  
  if (term == "nodeocov.log_gdp") {
    return("GDP per capita (outgoing)")
  }
  
  if (term == "nodeicov.log_gdp") {
    return("GDP per capita (incoming)")
  }
  
  if (term == "nodeocov.log_labor") {
    return("Labor force (outgoing)")
  }
  
  if (term == "nodeicov.log_labor") {
    return("Labor force (incoming)")
  }
  
  if (term == "edgecov.log_dist_144") {
    return("Log distance")
  }
  
  if (term == "edgecov.comlang_144") {
    return("Common language")
  }
  
  if (term == "edgecov.colony_144") {
    return("Colonial relationship")
  }
  
  if (term == "edgecov.contig_144") {
    return("Common border")
  }

  return(term)
}


# Geweke diagnosztika feldolgozása
parse_geweke <- function(
    diag_lines,
    terms
) {
  
  individual_idx <- grep(
    "^Individual P-values",
    diag_lines
  )
  
  joint_idx <- grep(
    "^Joint P-value",
    diag_lines
  )
  
  if (length(individual_idx) == 0 ||
      length(joint_idx) == 0) {
    
    return(
      list(
        table = data.frame(),
        joint_p = NA_real_,
        sample_size = NA_real_,
        thinning = NA_real_
      )
    )
  }
  
  i1 <- individual_idx[1]
  i2 <- joint_idx[1]
  
  # Csak a p-értékeket tartalmazó sorok
  p_lines <- diag_lines[
    (i1 + 1):(i2 - 1)
  ]
  
  p_lines <- p_lines[
    grepl(
      "^\\s*[0-9]",
      p_lines
    )
  ]
  
  p_values <- unlist(
    lapply(
      p_lines,
      function(x) {
        
        suppressWarnings(
          as.numeric(
            strsplit(
              trimws(x),
              "\\s+"
            )[[1]]
          )
        )
      }
    )
  )
  
  p_values <- p_values[
    !is.na(p_values)
  ]
  
  n_use <- min(
    length(p_values),
    length(terms)
  )
  
  geweke_table <- data.frame(
    Term = terms[1:n_use],
    Geweke_p = p_values[1:n_use],
    stringsAsFactors = FALSE
  )
  
  # Joint p
  joint_line <- diag_lines[
    joint_idx[1]
  ]
  
  joint_p <- suppressWarnings(
    as.numeric(
      sub(
        ".*Joint P-value.*:\\s*([0-9.eE+-]+).*",
        "\\1",
        joint_line
      )
    )
  )
  
  # sample size
  ss_line <- grep(
    "Sample size per chain =",
    diag_lines,
    value = TRUE
  )
  
  if (length(ss_line) == 0) {
    sample_size <- NA_real_
  } else {
    sample_size <- suppressWarnings(
      as.numeric(
        sub(
          ".*Sample size per chain =\\s*([0-9.eE+-]+).*",
          "\\1",
          ss_line[1]
        )
      )
    )
  }
  
  # thinning
  thinning_line <- grep(
    "Thinning interval =",
    diag_lines,
    value = TRUE
  )
  
  if (length(thinning_line) == 0) {
    thinning <- NA_real_
  } else {
    thinning <- suppressWarnings(
      as.numeric(
        sub(
          ".*Thinning interval =\\s*([0-9.eE+-]+).*",
          "\\1",
          thinning_line[1]
        )
      )
    )
  }
  
  list(
    table = geweke_table,
    joint_p = joint_p,
    sample_size = sample_size,
    thinning = thinning
  )
}


# GOF blokk feldolgozása
parse_gof_section <- function(
    gof_lines,
    section_title,
    section_type,
    year
) {
  
  idx <- grep(
    section_title,
    gof_lines,
    fixed = TRUE
  )
  
  if (length(idx) == 0) {
    return(data.frame())
  }
  
  start <- idx[1] + 1
  
  next_sections <- grep(
    "^Goodness-of-fit for ",
    gof_lines
  )
  
  end_candidates <- next_sections[
    next_sections > idx[1]
  ]
  
  if (length(end_candidates) == 0) {
    end <- length(gof_lines)
  } else {
    end <- min(end_candidates) - 1
  }
  
  block <- gof_lines[
    start:end
  ]
  
  # Sorok megtartása, ahol a 2. token numerikus
  parsed <- lapply(
    block,
    function(line) {
      
      parts <- strsplit(
        trimws(line),
        "\\s+"
      )[[1]]
      
      if (length(parts) < 6) {
        return(NULL)
      }
      
      observed <- suppressWarnings(
        as.numeric(parts[2])
      )
      
      if (is.na(observed)) {
        return(NULL)
      }
      
      data.frame(
        Year = year,
        Section = section_type,
        Statistic = parts[1],
        Observed = as.numeric(parts[2]),
        Sim_Min = as.numeric(parts[3]),
        Sim_Mean = as.numeric(parts[4]),
        Sim_Max = as.numeric(parts[5]),
        MC_p_value = as.numeric(parts[6]),
        stringsAsFactors = FALSE
      )
    }
  )
  
  bind_rows(parsed)
}


###############################################################################
# 5. STATIKUS 144x144 DYADICUS MÁTRIXOK
###############################################################################

comlang_144 <- make_block144(
  comlang,
  node_order,
  n_layers
)

log_dist_144 <- make_block144(
  log_dist,
  node_order,
  n_layers
)

colony_144 <- make_block144(
  colony,
  node_order,
  n_layers
)

contig_144 <- make_block144(
  contig,
  node_order,
  n_layers
)

###############################################################################
# 6. EREDMÉNYLISTÁK
###############################################################################

results_coefficients <- list()
results_model_fit <- list()
results_mcmc <- list()
results_geweke <- list()

results_gof_model <- list()
results_gof_degree <- list()
results_gof_triads <- list()
results_gof_layers <- list()

results_layers <- list()

results_checks <- list()

###############################################################################
# CONSTRAINED NESTED LAYER-EK
###############################################################################

# Ugyanaz a nested sample-space constraint, mint a korábban validált modellben.
nested_constraint <- ~ChangeStats(
  fix =
    ~L(~edges, ~`2` & !`1`) +
     L(~edges, ~`3` & !`2`) +
     L(~edges, ~`4` & !`3`)
)

###############################################################################
# FIT HELPER - constraint támogatással
###############################################################################

# Az eredeti fit_ergm_logged változata, egyetlen szükséges módosítással:
# az ergm() hívás megkapja a nested constraintet.
fit_ergm_logged <- function(formula, log_file, constraints = NULL) {
  
  message_file <- sub(
    "\\.log$",
    "_messages.log",
    log_file
  )
  
  warnings <- character()
  err <- NULL
  
  # A naplófájlok létrehozása és fejlécük mentése.
  writeLines(
    c(
      "ERGM fit log",
      paste("Started:", Sys.time()),
      ""
    ),
    con = log_file
  )
  
  writeLines(
    c(
      "ERGM message log",
      paste("Started:", Sys.time()),
      ""
    ),
    con = message_file
  )
  
  # Normál output: a Console-ban is látható marad.
  sink(log_file, append = TRUE, split = TRUE)
  
  # Üzenetek külön fájlba.
  msg_con <- file(message_file, open = "at")
  sink(msg_con, type = "message")
  
  model <- tryCatch(
    
    withCallingHandlers({
      
      cat("\nStarting ergm()...\n")
      
      ergm(
        formula,
        constraints = constraints,
        control = control.ergm(
          MCMC.burnin = 32768,
          MCMC.interval = 1024,
          MCMC.samplesize = 5000
        )
      )
      
    },
    
    warning = function(w) {
      
      warnings <<- c(
        warnings,
        conditionMessage(w)
      )
      
      cat(
        "\nWARNING:",
        conditionMessage(w),
        "\n"
      )
      
      invokeRestart("muffleWarning")
    }),
    
    error = function(e) {
      
      err <<- conditionMessage(e)
      
      cat(
        "\nERROR:",
        err,
        "\n"
      )
      
      NULL
    },
    
    finally = {
      
      # Először az üzenet-, utána a normál output-sinket zárjuk.
      sink(type = "message")
      close(msg_con)
      sink()
    }
  )
  
  # Mindkét napló beolvasása a parser számára.
  fit_output <- c(
    readLines(log_file, warn = FALSE),
    readLines(message_file, warn = FALSE)
  )
  
  if (!is.null(err)) {
    stop(err)
  }
  
  list(
    model = model,
    output = fit_output,
    warnings = unique(warnings),
    log_file = log_file,
    message_file = message_file
  )
}

###############################################################################
# 7. ÉVES CIKLUS
###############################################################################

for (year in years) {
  
  year_start <- Sys.time()
  
  cat("\n\n")
  cat("=================================================================\n")
  cat("YEAR:", year, "\n")
  cat("=================================================================\n")
  
  tic(paste("ERGM year", year))
  
  year_status <- "OK"
  year_error <- NA_character_
  
  ###########################################################################
  # 7.1 TRADE -> LAYEREK
  ###########################################################################
  
  trade_year <- trade[[as.character(year)]]
  
  # Trade check
  trade_ok <- ( # akkor true ha mind a 4 feltétel teljesül
    identical(rownames(trade_year), node_order) &&
      identical(colnames(trade_year), node_order) &&
      !anyNA(trade_year) &&
      all(trade_year >= 0)
  )
  
  if (!trade_ok) {
    stop(
      paste0(
        "Trade check failed in year ",
        year
      )
    )
  }
  
  trade_year_log <- log(
    trade_year + 1
  )
  
  positive_trade <- trade_year_log[
    trade_year_log > 0
  ]
  
  stopifnot(
    length(positive_trade) > 0
  )
  
  # ÉVES RELATÍV KÜSZÖBÖK
  # Ugyanaz a 50 / 75 / 90%-os felosztás, mint a legutóbb
  # validált constrained nested 2015-ös modellben.
  q50 <- as.numeric(
    quantile(
      positive_trade,
      probs = .50
    )
  )
  
  q75 <- as.numeric(
    quantile(
      positive_trade,
      probs = .75
    )
  )
  
  q90 <- as.numeric(
    quantile(
      positive_trade,
      probs = .90
    )
  )
  
  # Constrained nested layers
  L1 <- ifelse(
    trade_year_log > 0,
    1, 0
  )
  
  L2 <- ifelse(
    trade_year_log >= q50,
    1, 0
  )
  
  L3 <- ifelse(
    trade_year_log >= q75,
    1, 0
  )
  
  L4 <- ifelse(
    trade_year_log >= q90,
    1, 0
  )
  
  ###########################################################################
  # 7.2 LAYER CHECKS
  ###########################################################################
  
  # Nested rétegekben az átfedés szándékos.
  layer_overlap <- sum(
    L1 + L2 + L3 + L4 > 1
  )
  
  positive_dyads <- sum(
    trade_year > 0
  )
  
  layer_edge_counts <- c(
    L1_edges = sum(L1),
    L2_edges = sum(L2),
    L3_edges = sum(L3),
    L4_edges = sum(L4)
  )
  
  # Az összes pozitív kapcsolatnak L1-ben szerepelnie kell.
  layer_edge_sum_ok <- (
    sum(L1) == positive_dyads
  )
  
  # Nested hierarchia ellenőrzése: L4 ⊆ L3 ⊆ L2 ⊆ L1
  nested_layers_ok <- (
    all(L2 <= L1) &&
      all(L3 <= L2) &&
      all(L4 <= L3)
  )
  
  layer_mutual_counts <- c(
    L1_mutual_dyads = mutual_dyads(L1),
    L2_mutual_dyads = mutual_dyads(L2),
    L3_mutual_dyads = mutual_dyads(L3),
    L4_mutual_dyads = mutual_dyads(L4)
  )
  
  stopifnot(
    nested_layers_ok,
    layer_edge_sum_ok
  )
  
  ###########################################################################
  # 7.3 LAYER MATRICES MENTÉSE
  ###########################################################################
  
  saveRDS(
    list(
      Year = year,
      L1 = L1,
      L2 = L2,
      L3 = L3,
      L4 = L4,
      q50 = q50,
      q75 = q75,
      q90 = q90
    ),
    file = file.path(
      "annual_layers",
      paste0(
        "layers_",
        year,
        ".rds"
      )
    )
  )
  
  # összefoglaló tábla minden évre - a végén egy 10 soros összefoglaló tábla mert egymás alá fűzi
  results_layers[[length(results_layers) + 1]] <- data.frame(
    Year = year,
    q50 = q50,
    q75 = q75,
    q90 = q90,
    L1_edges = sum(L1),
    L2_edges = sum(L2),
    L3_edges = sum(L3),
    L4_edges = sum(L4),
    L1_mutual_dyads = mutual_dyads(L1),
    L2_mutual_dyads = mutual_dyads(L2),
    L3_mutual_dyads = mutual_dyads(L3),
    L4_mutual_dyads = mutual_dyads(L4),
    Positive_dyads = positive_dyads,
    Layer_overlap = layer_overlap
  )
  
  ###########################################################################
  # 7.5 NETWORK OBJECTUMOK
  ###########################################################################
  
  L1_net <- network(
    L1,
    directed = TRUE,
    matrix.type = "adjacency"
  )
  
  L2_net <- network(
    L2,
    directed = TRUE,
    matrix.type = "adjacency"
  )
  
  L3_net <- network(
    L3,
    directed = TRUE,
    matrix.type = "adjacency"
  )
  
  L4_net <- network(
    L4,
    directed = TRUE,
    matrix.type = "adjacency"
  )
  
  # Node-sorrend ellenőrzése
  stopifnot(
    identical(
      network.vertex.names(L1_net),
      node_order
    ),
    identical(
      network.vertex.names(L2_net),
      node_order
    ),
    identical(
      network.vertex.names(L3_net),
      node_order
    ),
    identical(
      network.vertex.names(L4_net),
      node_order
    )
  )
  
  LL_trade <- Layer(
    L1 = L1_net,
    L2 = L2_net,
    L3 = L3_net,
    L4 = L4_net
  )
  
  # ellenőrzés: valóban azt a 4 megfigyelt mátrixot kapom vissza amit a modellbe betettem
  observed_block_check <- extract_layer_matrices(LL_trade)
  
  stopifnot(
    isTRUE(all.equal(
      unname(observed_block_check$L1),
      unname(L1)
    )),
    isTRUE(all.equal(
      unname(observed_block_check$L2),
      unname(L2)
    )),
    isTRUE(all.equal(
      unname(observed_block_check$L3),
      unname(L3)
    )),
    isTRUE(all.equal(
      unname(observed_block_check$L4),
      unname(L4)
    ))
  )
  
  rm(observed_block_check)
  ###########################################################################
  # 7.6 WORLD BANK NODE COVARIATES
  ###########################################################################
  
  wb_year <- wb[[as.character(year)]]
  
  log_hightech_144 <- rep(
    wb_year$log_hightech,
    4
  )
  
  log_gdp_144 <- rep(
    wb_year$log_gdp,
    4
  )
  
  log_labor_144 <- rep(
    wb_year$log_labor,
    4
  )
  
  LL_trade %v% "log_hightech" <- log_hightech_144
  LL_trade %v% "log_gdp" <- log_gdp_144
  LL_trade %v% "log_labor" <- log_labor_144
  
  # Node covariate check
  same_across_layers <- function(x) {
    
    all(x[1:36] == x[37:72]) &&
      all(x[1:36] == x[73:108]) &&
      all(x[1:36] == x[109:144])
  }
  
  node_covariates_ok <- (
    same_across_layers(
      LL_trade %v% "log_hightech"
    ) &&
      same_across_layers(
        LL_trade %v% "log_gdp"
      ) &&
      same_across_layers(
        LL_trade %v% "log_labor"
      )
  )
  
  stopifnot(node_covariates_ok)
  
  ###########################################################################
  # 7.7 MODEL 8 FORMULA
  ###########################################################################
  
  model_formula <- LL_trade ~
    
    # Layer-specific edge density
    L(~edges, ~`1`) +
    L(~edges, ~`2`) +
    L(~edges, ~`3`) +
    L(~edges, ~`4`) +
    
    # Layer-specific mutuality
    mutualL(
      Ls = list(
        ~`1`,
        ~`1`
      )
    ) +
    
    mutualL(
      Ls = list(
        ~`2`,
        ~`2`
      )
    ) +
    
    mutualL(
      Ls = list(
        ~`3`,
        ~`3`
      )
    ) +
    
    mutualL(
      Ls = list(
        ~`4`,
        ~`4`
      )
    ) +
    
    # Node-level covariates
    nodeocov("log_hightech") +
    nodeicov("log_hightech") +
    
    nodeocov("log_gdp") +
    nodeicov("log_gdp") +
    
    nodeocov("log_labor") +
    nodeicov("log_labor") +
    
    # Dyadic controls
    edgecov(log_dist_144) +
    edgecov(comlang_144) +
    edgecov(colony_144) +
    edgecov(contig_144) 
  
  ###########################################################################
  # 7.8 MODEL FIT
  ###########################################################################
  
  # minden évnek saját reprodukálható seedje lesz (12015, ... , 12024)
  set.seed(
    10000 + year
  )
  
  log_file <- file.path(
    "annual_logs",
    paste0(
      "ergm_ml8_",
      year,
      ".log"
    )
  )
  
  # a tényleges ergm-futtatás
  fit <- fit_ergm_logged(
    model_formula,
    log_file,
    constraints = nested_constraint
  )
  
  model_ml_8 <- fit$model
  
  mcmle_info <- parse_mcmle_info(
    fit$output
  )
  
  fit_metrics <- parse_model_fit(
    model_ml_8
  )
  
  ###########################################################################
  # 7.9 COEFFICIENTEK
  ###########################################################################
  
  # a frissen lefuttatott 1 év modelljének együtthatóinak elrendezése és eltárolása
  coef_year <- extract_coefficients(
    model_ml_8,
    year
  )
  
  # labelek gyártása hogy olvashatóbb legyen
  coef_year$Group <- vapply(
    coef_year$Term,
    classify_term,
    character(1)
  )
  
  coef_year$Label <- vapply(
    coef_year$Term,
    term_label,
    character(1)
  )
  
  results_coefficients[[length(results_coefficients) + 1]] <-
    coef_year
  
  ###########################################################################
  # 7.10 MODEL FIT TÁBLA
  ###########################################################################
  
  results_model_fit[[length(results_model_fit) + 1]] <-
    data.frame(
      Year = year,
      AIC = fit_metrics$AIC,
      BIC = fit_metrics$BIC,
      Null_Deviance = fit_metrics$Null_Deviance,
      Null_df = fit_metrics$Null_df,
      Residual_Deviance = fit_metrics$Residual_Deviance,
      Residual_df = fit_metrics$Residual_df,
      MC_SE = fit_metrics$MC_SE,
      MCMLE_convergence_p =
        mcmle_info$MCMLE_convergence_p,
      MCMLE_convergence_p_reported =
        mcmle_info$MCMLE_convergence_p_reported,
      MCMLE_converged_99 =
        mcmle_info$MCMLE_converged_99,
      MCMLE_iterations =
        mcmle_info$MCMLE_iterations,
      Warning_count =
        length(fit$warnings),
      Warnings =
        ifelse(
          length(fit$warnings) == 0,
          NA_character_,
          paste(
            fit$warnings,
            collapse = " | "
          )
        )
    )
  
  ###########################################################################
  # 7.11 MCMC DIAGNOSTICS
  ###########################################################################
  
  mcmc_pdf <- file.path(
    "model diagnostics",
    paste0(
      "mcmc_diagnostics_ml8_",
      year,
      ".pdf"
    )
  )
  
  pdf(
    mcmc_pdf,
    width = 10,
    height = 8
  )
  
  diag_lines <- capture.output(
    mcmc.diagnostics(
      model_ml_8
    )
  )
  
  dev.off()
  
  geweke_info <- parse_geweke(
    diag_lines,
    coef_year$Term
  )
  
  geweke_year <- geweke_info$table
  
  if (nrow(geweke_year) > 0) {
    
    geweke_year$Year <- year
    
    geweke_year$Group <- vapply(
      geweke_year$Term,
      classify_term,
      character(1)
    )
    
    geweke_year$Label <- vapply(
      geweke_year$Term,
      term_label,
      character(1)
    )
  }
  
  results_geweke[[length(results_geweke) + 1]] <-
    geweke_year
  
  if (nrow(geweke_year) > 0) {
    
    min_idx <- which.min(
      geweke_year$Geweke_p
    )
    
    min_geweke_p <-
      geweke_year$Geweke_p[min_idx]
    
    min_geweke_term <-
      geweke_year$Term[min_idx]
    
    n_geweke_bad <-
      sum(
        geweke_year$Geweke_p < 0.05,
        na.rm = TRUE
      )
    
  } else {
    
    min_geweke_p <- NA_real_
    min_geweke_term <- NA_character_
    n_geweke_bad <- NA_integer_
  }
  
  # rövid éves összefoglaló MCMC-diagnosztikai tábla
  results_mcmc[[length(results_mcmc) + 1]] <-
    data.frame(
      Year = year,
      Joint_Geweke_p =
        geweke_info$joint_p,
      Min_Geweke_p =
        min_geweke_p,
      Min_Geweke_term =
        min_geweke_term,
      N_Geweke_p_below_0.05 =
        n_geweke_bad,
      MCMC_sample_size =
        geweke_info$sample_size,
      Thinning_interval =
        geweke_info$thinning
    )
  
  ###########################################################################
  # 7.12 GOF
  ###########################################################################
  
  set.seed(
    20000 + year
  )
  
  gof_ml_8_full <- gof(
    model_ml_8,
    GOF = ~
      model +
      idegree +
      odegree +
      triadcensus
  )
  
  gof_pdf <- file.path(
    "annual_gof",
    paste0(
      "gof_ml8_",
      year,
      ".rds"
    )
  )
  
  saveRDS(
    gof_ml_8_full,
    gof_pdf
  )
  
  gof_lines <- capture.output(
    print(gof_ml_8_full)
  )
  
  # Model GOF
  gof_model_year <- parse_gof_section(
    gof_lines,
    "Goodness-of-fit for model statistics",
    "model",
    year
  )
  
  # Out-degree
  gof_outdegree_year <- parse_gof_section(
    gof_lines,
    "Goodness-of-fit for out-degree",
    "out-degree",
    year
  )
  
  # In-degree
  gof_indegree_year <- parse_gof_section(
    gof_lines,
    "Goodness-of-fit for in-degree",
    "in-degree",
    year
  )
  
  # Triad census
  gof_triad_year <- parse_gof_section(
    gof_lines,
    "Goodness-of-fit for triad census",
    "triadcensus",
    year
  )
  
  results_gof_model[[length(results_gof_model) + 1]] <-
    gof_model_year
  
  results_gof_degree[[length(results_gof_degree) + 1]] <-
    bind_rows(
      gof_outdegree_year,
      gof_indegree_year
    )
  
  results_gof_triads[[length(results_gof_triads) + 1]] <-
    gof_triad_year
  
  ###############################################################################
  # 7.12B RÉTEGENKÉNTI GOF SZIMULÁCIÓVAL
  ###############################################################################
  
  set.seed(30000 + year)
  
  layer_sim_networks <- simulate(
    model_ml_8,
    nsim = n_layer_gof_sim,
    output = "network"
  )
  
  layer_gof_year <- compute_layer_gof(
    observed_layers = list(
      L1 = L1,
      L2 = L2,
      L3 = L3,
      L4 = L4
    ),
    sim_networks = layer_sim_networks,
    year = year
  )
  
  # A részletes, rétegenkénti GOF mentése az adott évre
  results_gof_layers[[length(results_gof_layers) + 1]] <-
    layer_gof_year
  
  saveRDS(
    layer_gof_year,
    file.path(
      "annual_gof",
      paste0("gof_ml8_layers_", year, ".rds")
    )
  )
  
  # Szükségtelen objektumok felszabadítása
  rm(layer_sim_networks, layer_gof_year)
  
  ###########################################################################
  # 7.13 ROUTE GOF / DEGREE / TRIAD DIAGNOSZTIKAI ÖSSZEFOGLALÓ
  ###########################################################################
  
  if (nrow(gof_model_year) > 0) {
    n_gof_model_bad <-
      sum(
        gof_model_year$MC_p_value < 0.05,
        na.rm = TRUE
      )
  } else {
    n_gof_model_bad <- NA_integer_
  }
  
  degree_year <- bind_rows(
    gof_outdegree_year,
    gof_indegree_year
  )
  
  if (nrow(degree_year) > 0) {
    n_gof_degree_bad <-
      sum(
        degree_year$MC_p_value < 0.05,
        na.rm = TRUE
      )
  } else {
    n_gof_degree_bad <- NA_integer_
  }
  
  if (nrow(gof_triad_year) > 0) {
    n_gof_triad_bad <-
      sum(
        gof_triad_year$MC_p_value < 0.05,
        na.rm = TRUE
      )
  } else {
    n_gof_triad_bad <- NA_integer_
  }
  
  ###########################################################################
  # 7.14 CHECKS
  ###########################################################################
  
  results_checks[[length(results_checks) + 1]] <-
    data.frame(
      Year = year,
      
      Trade_ok = trade_ok,
      Node_covariates_ok =
        node_covariates_ok,
      
      Layer_overlap = layer_overlap,
      Nested_layers_ok =
        nested_layers_ok,
      Layer_edge_sum_ok =
        layer_edge_sum_ok,
      
      CEPII_no_NA =
        !anyNA(comlang) &&
        !anyNA(log_dist) &&
        !anyNA(colony) &&
        !anyNA(contig),
      
      MCMLE_converged_99 =
        mcmle_info$MCMLE_converged_99,
      
      Joint_Geweke_ok =
        geweke_info$joint_p >= 0.05,
      
      Min_Geweke_p =
        min_geweke_p,
      
      N_Geweke_p_below_0.05 =
        n_geweke_bad,
      
      N_GOF_model_p_below_0.05 =
        n_gof_model_bad,
      
      N_GOF_degree_p_below_0.05 =
        n_gof_degree_bad,
      
      N_GOF_triad_p_below_0.05 =
        n_gof_triad_bad,
      
      Model_status = year_status,
      Error_message = year_error,
      
      stringsAsFactors = FALSE
    )
  
  ###########################################################################
  # 7.15 MODELL ÉS GOF MENTÉSE
  ###########################################################################
  
  saveRDS(
    model_ml_8,
    file.path(
      "annual_models",
      paste0(
        "model_ml8_",
        year,
        ".rds"
      )
    )
  )
  
  ###########################################################################
  # 7.16 ÉVES EREDMÉNYEK IDEIGLENES MENTÉSE
  ###########################################################################
  
  saveRDS(
    list(
      coefficients =
        bind_rows(results_coefficients),
      model_fit =
        bind_rows(results_model_fit),
      mcmc =
        bind_rows(results_mcmc),
      geweke =
        bind_rows(results_geweke),
      gof_model =
        bind_rows(results_gof_model),
      gof_degree =
        bind_rows(results_gof_degree),
      gof_triads =
        bind_rows(results_gof_triads),
      gof_layers =
        bind_rows(results_gof_layers),
      layers =
        bind_rows(results_layers),
      checks =
        bind_rows(results_checks)
    ),
    file.path(
      "annual_results",
      paste0(
        "results_up_to_",
        year,
        ".rds"
      )
    )
  )
  
  ###########################################################################
  # 7.17 ÉVES FUTÁSI IDŐ
  ###########################################################################
  
  runtime_seconds <- as.numeric(
    difftime(
      Sys.time(),
      year_start,
      units = "secs"
    )
  )
  
  # Futási időt hozzáadjuk a model fit táblához
  results_model_fit[[length(results_model_fit)]]$Runtime_seconds <-
    runtime_seconds
  
  # Ellenőrzés kedvéért frissítjük az utolsó modell-fit sort
  cat(
    "\nYear",
    year,
    "finished in",
    round(runtime_seconds, 2),
    "seconds.\n"
  )
  
  cat(
    "AIC:",
    fit_metrics$AIC,
    "| BIC:",
    fit_metrics$BIC,
    "| Joint Geweke:",
    geweke_info$joint_p,
    "\n"
  )
  
  if (length(fit$warnings) > 0) {
    
    cat(
      "Warnings:",
      paste(
        fit$warnings,
        collapse = " | "
      ),
      "\n"
    )
  }
  
  # Memória felszabadítása
  rm(
    model_ml_8,
    gof_ml_8_full,
    LL_trade,
    L1_net,
    L2_net,
    L3_net,
    L4_net
  )
  
  gc()
  
  toc()
}


###############################################################################
# 8. ÖSSZES EREDMÉNY ÖSSZEÁLLÍTÁSA
###############################################################################

# listákba gyűjti évenként az eredményeket
results_coefficients_df <- bind_rows(
  results_coefficients
)

results_model_fit_df <- bind_rows(
  results_model_fit
)

results_mcmc_df <- bind_rows(
  results_mcmc
)

results_geweke_df <- bind_rows(
  results_geweke
)

results_gof_model_df <- bind_rows(
  results_gof_model
)

results_gof_degree_df <- bind_rows(
  results_gof_degree
)

results_gof_triads_df <- bind_rows(
  results_gof_triads
)

results_layers_df <- bind_rows(
  results_layers
)

results_checks_df <- bind_rows(
  results_checks
)

results_gof_layers_df <- bind_rows(
  results_gof_layers
)

results_gof_layer_summary_df <- results_gof_layers_df %>%
  group_by(
    Year,
    Layer,
    Statistic
  ) %>%
  summarise(
    N_categories = n(),
    N_MC_p_below_0_05 = sum(
      MC_two_sided_p < 0.05,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

###############################################################################
# 10. ÉVES ÁTTEKINTŐ TÁBLA
###############################################################################

annual_overview_df <- results_model_fit_df %>%
  select(
    Year,
    AIC,
    BIC,
    Null_Deviance,
    Residual_Deviance,
    MC_SE,
    MCMLE_iterations,
    MCMLE_convergence_p_reported,
    MCMLE_converged_99,
    Runtime_seconds
  ) %>%
  left_join(
    results_mcmc_df %>%
      select(
        Year,
        Joint_Geweke_p,
        Min_Geweke_p,
        Min_Geweke_term,
        N_Geweke_p_below_0.05,
        MCMC_sample_size,
        Thinning_interval
      ),
    by = "Year"
  ) %>%
  left_join(
    results_layers_df %>%
      select(
        Year,
        q50,
        q75,
        q90,
        L1_edges,
        L2_edges,
        L3_edges,
        L4_edges,
        L1_mutual_dyads,
        L2_mutual_dyads,
        L3_mutual_dyads,
        L4_mutual_dyads,
        Positive_dyads
      ),
    by = "Year"
  )

###############################################################################
# 12. FINAL ELLENŐRZŐ KIÍRÁS
###############################################################################

cat("\n\n")
cat("=================================================================\n")
cat("ANNUAL NESTED MODEL 8 RUN FINISHED\n")
cat("=================================================================\n")

cat(
  "Years successfully processed:",
  paste(
    unique(
      results_model_fit_df$Year
    ),
    collapse = ", "
  ),
  "\n"
)

cat(
  "Number of coefficient rows:",
  nrow(results_coefficients_df),
  "\n"
)

cat(
  "Number of GOF model rows:",
  nrow(results_gof_model_df),
  "\n"
)

###############################################################################
# 13. EXCEL EREDMÉNYFÁJL
###############################################################################

output_xlsx <- "ml_ergm_nested_model8_results.xlsx"

writexl::write_xlsx(
  list(
    annual_overview = annual_overview_df,
    coefficients = results_coefficients_df,
    model_fit = results_model_fit_df,
    mcmc = results_mcmc_df,
    geweke = results_geweke_df,
    gof_model = results_gof_model_df,
    gof_degree = results_gof_degree_df,
    gof_triads = results_gof_triads_df,
    gof_layers = results_gof_layers_df,
    gof_layer_summary = results_gof_layer_summary_df,
    layers = results_layers_df,
    checks = results_checks_df
  ),
  path = output_xlsx
)

cat(
  "\nExcel eredményfájl létrehozva:",
  output_xlsx,
  "\n"
)

###############################################################################
# VÉGE
###############################################################################

getwd()

file.info(
  file.path("annual_logs", "ergm_ml8_2024.log")
)[, c("size", "mtime")]

readLines(
  file.path("annual_logs", "ergm_ml8_2024.log"),
  n = 40,
  warn = FALSE
)

readLines(
  file.path("annual_logs", "ergm_ml8_2024_messages.log"),
  n = 40,
  warn = FALSE
)


