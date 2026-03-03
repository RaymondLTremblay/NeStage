# R/nestage-logging.R

#' Append a one-line entry to logs/CHANGELOG.md
#'
#' @param message Text to log.
#' @param log_path Path to changelog (default 'logs/CHANGELOG.md').
#' @return Invisibly returns the appended line.
#' @export
nestage_log_run <- function(message, log_path = "logs/CHANGELOG.md") {
  ver <- nestage_version()
  stamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
  line <- sprintf("%s | %s | %s\n", stamp, ver, message)
  dir.create(dirname(log_path), showWarnings = FALSE, recursive = TRUE)
  cat(line, file = log_path, append = TRUE)
  invisible(line)
}
