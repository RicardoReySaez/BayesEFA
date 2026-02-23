#' @importFrom utils packageVersion
#' @noRd

.onAttach <- function(libname, pkgname) {
  # ───────────────────────── #
  #    Get package version    #
  # ───────────────────────── #

  pkg_version <- tryCatch(
    as.character(utils::packageVersion("BayesEFA")),
    error = function(e) "unknown"
  )

  # ─────────────────────────── #
  #    Build startup message    #
  # ─────────────────────────── #

  line_width <- 60
  thick_line <- paste0(rep("\u2500", line_width), collapse = "")

  header <- sprintf(" BayesEFA (v%s): Bayesian Exploratory Factor Analysis", pkg_version)

  # Check for cmdstanr
  has_cmdstanr <- requireNamespace("cmdstanr", quietly = TRUE)

  msgs <- c()

  if (!has_cmdstanr) {
    # If cmdstanr is missing, recommend it
    msgs <- c(
      msgs,
      " * Recommendation: We strongly recommend installing 'cmdstanr'.",
      "   It provides the latest Stan version and faster estimation.",
      "   Install instructions: https://mc-stan.org/cmdstanr/articles/cmdstanr.html"
    )
  } else {
    # If cmdstanr is present, check for cmdstan
    options(cmdstanr_warn_inits = FALSE)
    cmdstan_path <- tryCatch(cmdstanr::cmdstan_path(), error = function(e) NULL)

    if (is.null(cmdstan_path)) {
      # cmdstanr installed but no cmdstan found
      msgs <- c(
        msgs,
        " * Note: 'cmdstanr' (the R package) is installed but no CmdStan installation was found.",
        "   Please clean/install CmdStan to enable the faster backend.",
        "   See 'Installing CmdStan' at: https://mc-stan.org/cmdstanr/articles/cmdstanr.html"
      )
    }
  }

  # Assemble message
  msg_parts <- c(
    "",
    thick_line,
    header,
    if (length(msgs) > 0) {
      c(paste0(rep("\u2500", line_width), collapse = ""), msgs)
    } else {
      NULL
    },
    thick_line,
    ""
  )

  packageStartupMessage(paste(msg_parts, collapse = "\n"))
}
