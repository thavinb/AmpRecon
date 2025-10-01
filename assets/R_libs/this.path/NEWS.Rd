\name{NEWS}
%% used in ./inst/NEWS.Rd
\newcommand{\codeorsQuote}{\ifelse{text}{\sQuote{#1}}{\code{#1}}}
\newcommand{\commandorsQuote}{\ifelse{text}{\sQuote{#1}}{\command{#1}}}
\newcommand{\optionorsQuote}{\ifelse{text}{\sQuote{#1}}{\option{#1}}}
\newcommand{\pkgorsQuote}{\ifelse{text}{\sQuote{#1}}{\pkg{#1}}}


%% does not work well with characters that have
%% special meaning in HTML or LaTeX or Rd
\newcommand{\string}{\ifelse{text}{"#1"}{\code{"#1"}}}
\newcommand{\blankstring}{\ifelse{text}{""}{\code{""}}}
\newcommand{\codelink2}{\ifelse{html}{\code{\link[#1]{#2}}}{\codeorsQuote{#2}}}
\newcommand{\codelink3}{\ifelse{html}{\code{\link[#1]{#2}#3}}{\codeorsQuote{#2#3}}}
\newcommand{\codelink4}{\ifelse{html}{\code{#1\link[#2]{#3}#4}}{\codeorsQuote{#1#3#4}}}


\newcommand{\CRANpackage}{\href{https://CRAN.R-project.org/package=#1}{\pkgorsQuote{package:#1}}}


\newcommand{\Emacs}{\href{https://www.gnu.org/software/emacs/}{Emacs}}
\newcommand{\ESS}{\href{https://ess.r-project.org/}{ESS}}
\newcommand{\Jupyter}{\href{https://jupyter.org/}{Jupyter}}
\newcommand{\Positron}{\href{https://github.com/posit-dev/positron}{Positron}}
\newcommand{\Python}{\href{https://www.python.org/}{Python}}
\newcommand{\radian}{\href{https://github.com/randy3k/radian}{radian}}
\newcommand{\REditorSupport}{\href{https://code.visualstudio.com/docs/languages/r}{REditorSupport}}
\newcommand{\RStudio}{\href{https://posit.co/products/open-source/rstudio/}{RStudio}}
\newcommand{\VSCode}{\href{https://code.visualstudio.com/}{VSCode}}


\newcommand{\sourcelink}{\codelink3{base}{source}{()}}
\newcommand{\syssourcelink}{\codelink3{base}{sys.source}{()}}
\newcommand{\debugSourceinRStudiolink}{\ifelse{text}{\sQuote{debugSource()} in \sQuote{RStudio}}{\href{https://support.posit.co/hc/en-us/articles/205612627-Debugging-with-the-RStudio-IDE}{\code{debugSource}}\code{()} in \sQuote{\RStudio}}}
\newcommand{\compilerloadcmplink}{\codelink3{compiler:compile}{compiler::loadcmp}{()}}
\newcommand{\utilsSweavelink}{\codelink3{utils:Sweave}{utils::Sweave}{()}}
\newcommand{\boxuselink}{\codelink3{box:use}{box::use}{()}}
\newcommand{\knitrknitlink}{\codelink3{knitr:knit}{knitr::knit}{()}}
\newcommand{\plumberplumblink}{\codelink3{plumber:plumb}{plumber::plumb}{()}}
\newcommand{\shinyrunApplink}{\codelink3{shiny:runApp}{shiny::runApp}{()}}
\newcommand{\targetslink}{\CRANpackage{targets}}
\newcommand{\targetstarcallrinnertrylink}{\codelink3{targets:tar_callr_inner_try}{targets::tar_callr_inner_try}{()}}
\newcommand{\targetstarloadglobalslink}{\codelink3{targets:tar_load_globals}{targets::tar_load_globals}{()}}
\newcommand{\targetstarsourcelink}{\codelink3{targets:tar_source}{targets::tar_source}{()}}
\newcommand{\targetstarworkspacelink}{\codelink3{targets:tar_workspace}{targets::tar_workspace}{()}}
\newcommand{\testthatsourcefilelink}{\codelink3{testthat:source_file}{testthat::source_file}{()}}
\newcommand{\Rscriptlink}{\codelink2{utils}{Rscript}}

\newcommand{\sourcenolink}{\codeorsQuote{source()}}
\newcommand{\syssourcenolink}{\codeorsQuote{sys.source()}}
\newcommand{\debugSourceinRStudionolink}{\codeorsQuote{debugSource()} in \sQuote{RStudio}}
\newcommand{\compilerloadcmpnolink}{\codeorsQuote{compiler::loadcmp()}}
\newcommand{\utilsSweavenolink}{\codeorsQuote{utils::Sweave()}}
\newcommand{\boxusenolink}{\codeorsQuote{box::use()}}
\newcommand{\knitrknitnolink}{\codeorsQuote{knitr::knit()}}
\newcommand{\plumberplumbnolink}{\codeorsQuote{plumber::plumb()}}
\newcommand{\shinyrunAppnolink}{\codeorsQuote{shiny::runApp()}}
\newcommand{\targetsnolink}{\pkgorsQuote{package:targets}}
\newcommand{\targetstarcallrinnertrynolink}{\codeorsQuote{targets::tar_callr_inner_try()}}
\newcommand{\targetstarloadglobalsnolink}{\codeorsQuote{targets::tar_load_globals()}}
\newcommand{\targetstarsourcenolink}{\codeorsQuote{targets::tar_source()}}
\newcommand{\targetstarworkspacenolink}{\codeorsQuote{targets::tar_workspace()}}
\newcommand{\testthatsourcefilenolink}{\codeorsQuote{testthat::source_file()}}
\newcommand{\Rscriptnolink}{\codeorsQuote{Rscript}}


\newcommand{\syspathcompatible}{\itemize{\item{\source#1} \item{\syssource#1} \item{\debugSourceinRStudio#1} \item{\compilerloadcmp#1} \item{\utilsSweave#1} \item{\boxuse#1} \item{\knitrknit#1} \item{\plumberplumb#1} \item{\shinyrunApp#1} \item{\targets#1} \item{\testthatsourcefile#1}}}
\newcommand{\syspathcompatiblelinks}{\syspathcompatible{link}}
\newcommand{\syspathcompatiblenolinks}{\syspathcompatible{nolink}}


\newcommand{\sources}{\itemize{\item{\source#1} \item{\syssource#1} \item{\debugSourceinRStudio#1} \item{\compilerloadcmp#1} \item{\utilsSweave#1} \item{\boxuse#1} \item{\knitrknit#1} \item{\plumberplumb#1} \item{\shinyrunApp#1} \item{\targetstarcallrinnertry#1 \cr \targetstarloadglobals#1 \cr \targetstarsource#1 \cr \targetstarworkspace#1} \item{\testthatsourcefile#1}}}
\newcommand{\sourceslinks}{\sources{link}}
\newcommand{\sourcesnolinks}{\sources{nolink}}
\title{\pkg{this.path} News}
\section{CHANGES IN VERSION 2.7.0 (2025-04-10)}{
  \subsection{SIGNIFICANT USER-VISIBLE CHANGES}{
    \itemize{
      \item{Previously defunct functions \code{set.sys.path.jupyter()} and
        \code{fileArgs()} were removed entirely.}

      \item{\code{path.join()} now respects the declared encodings of its
        arguments and appropriately marks the encodings of its output, similar
        to \code{paste()}.}

      \item{\code{path.split()} and \code{path.unsplit()} now respect the
        declared encodings of their arguments and appropriately mark the
        encodings of their outputs.}
    }
  }
  \subsection{NEW FEATURES}{
    \itemize{
      \item{\code{set.sys.path()} now accepts argument \codeorsQuote{delayed}
        indicating whether the normalizing of the path should be delayed.}

      \item{\code{make_fix_funs()} now accepts argument \codeorsQuote{delayed}
        indicating whether the normalizing of \codeorsQuote{file} should be
        delayed.}
    }
  }
}
\section{CHANGES IN VERSION 2.6.0 (2024-12-18)}{
  \subsection{SIGNIFICANT USER-VISIBLE CHANGES}{
    \itemize{
      \item{Previously defunct functions \code{reset.this.proj()},
        \code{set.this.path.jupyter()}, \code{inside.source()},
        \code{set.this.path()}, and \code{unset.this.path()} were removed
        entirely.}
    }
  }
  \subsection{NEW FEATURES}{
    \itemize{
      \item{\code{this.path()} is now compatible with GUI \sQuote{\Positron}.}

      \item{\code{with_script_path()} now accepts argument \codeorsQuote{file}
        specifying the pathname of the file.}

      \item{Added \code{make_fix_funs()}, an alias of \code{path.functions()}.}
    }
  }
  \subsection{BUG FIXES}{
    \itemize{
      \item{\code{this.path(contents = TRUE)} in \sQuote{Rgui} in which the
        source document contains no text now returns \code{character(0)}
        instead of the previously incorrect \blankstring.}

      \item{\code{this.path()} in \sQuote{\RStudio} and \sQuote{\VSCode} will
        now declare a \string{UTF-8} file path encoding only on Windows.}

      \item{Improved detection of \sQuote{\RStudio}, no longer relying on
        environment variable \string{RSTUDIO_SESSION_PID} which may be unset.}
    }
  }
}
\section{CHANGES IN VERSION 2.5.0 (2024-06-29)}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item{\code{progArgs()} now accepts argument \codeorsQuote{ifnotfound}
        specifying what to do when no program arguments are found.}
    }
  }
  \subsection{BUG FIXES}{
    \itemize{
      \item{\code{this.path(original = TRUE)} used inside a \string{gzfile},
        \string{bzfile}, \string{xzfile}, or \string{fifo} connection will no
        longer return a \samp{file://} URL.}

      \item{\code{this.path()} once again works in \sQuote{Rgui} running in a
        Japanese or Turkish locale.}
    }
  }
}
\section{CHANGES IN VERSION 2.4.0 (2024-02-16)}{
  \subsection{SIGNIFICANT USER-VISIBLE CHANGES}{
    \itemize{
      \item{\code{with_init.file()} now accepts a braced expression, the
        sub-expressions of which are evaluated and auto-printed.}
    }
  }
  \subsection{NEW FEATURES}{
    \itemize{
      \item{\code{this.path()} is now compatible with function
        \code{utils::Sweave()}.}

      \item{\code{this.path()} is now compatible with background jobs in
        \sQuote{\RStudio}.}

      \item{Added \code{with_script_path()} which sets the current script in
        \pkgorsQuote{package:box} to \code{this.path()}, then evaluates its
        argument, most commonly a \pkgorsQuote{package:box} import statement.}

      \item{Added \code{make_fix_file()} which makes a function that constructs
        file paths, starting with the project root of \code{this.dir()}.}
    }
  }
  \subsection{DEPRECATED AND DEFUNCT}{
    \itemize{
      \item{\code{set.sys.path.jupyter()} is defunct and should be replaced
        with \code{set.jupyter.path()}.}

      \item{\code{fileArgs()} is defunct and should be replaced
        with \code{progArgs()}.}
    }
  }
  \subsection{BUG FIXES}{
    \itemize{
      \item{\code{this.path()} now works correctly in the site-wide startup
        profile file.}

      \item{\code{this.path()} now works correctly with source file aliases.}

      \item{Updated \code{shFILE()} to work better with older versions of \R.}

      \item{Improved \code{LINENO()}, specifically \code{sys.LINENO()}, when
        comparing paths to determine the appropriate source reference.}

      \item{\code{this.proj()} now throws an error that it does not work for
        URL pathnames.}

      \item{Removed the restriction that \code{reset.proj()} could only be
        called from a top level context.}

      \item{\code{withArgs()} will now be found anywhere on the stack as long
        as it and the \sourcenolink-like function are called in the same frame.

        e.g. \code{withArgs(try(source("./example.R")), 1, 2, 3)} now works as
        intended.}

      \item{Unified behaviour of \code{is.main()} across all GUIs.}

      \item{All promises are now braced, giving a more sensible call stack in
        case of warnings / / errors.}

      \item{Namespaces loaded with \code{requireNamespace()} are now loaded
        quietly.}
    }
  }
}
\section{CHANGES IN VERSION 2.3.1 (2023-12-10)}{
  \subsection{BUG FIXES}{
    \itemize{
      \item{\code{sys.srcref()} now works when passed as a non-forced argument
        to a function containing source references. By extension,
        \code{this.path()} and \code{LINENO()} now work as well.}
    }
  }
}
\section{CHANGES IN VERSION 2.3.0 (2023-12-02)}{
  \subsection{SIGNIFICANT USER-VISIBLE CHANGES}{
    \itemize{
      \item{Previously defunct functions \code{this.path2()},
        \code{this.dir2()}, \code{this.dir3()}, and \code{normalized.shFILE()}
        were removed entirely.}
    }
  }
  \subsection{NEW FEATURES}{
    \itemize{
      \item{\code{this.path()} is now compatible with GUI \sQuote{\Emacs} +
        package \sQuote{\ESS}.}
    }
  }
}
\section{CHANGES IN VERSION 2.2.0 (2023-10-29)}{
  \subsection{SIGNIFICANT USER-VISIBLE CHANGES}{
    \itemize{
      \item{Active bindings \codeorsQuote{FILE} and \codeorsQuote{LINE} have
        been converted to regular bindings \code{FILE()} and \code{LINE()}
        after the unfortunate realization that active binding functions are
        always evaluated in the global environment.}
    }
  }
  \subsection{NEW FEATURES}{
    \itemize{
      \item{\code{this.path()} works significantly better with source
        references.}

      \item{\code{this.path()} is now compatible with package
        \CRANpackage{targets}.}

      \item{Added \code{with_site.file()} and \code{with_init.file()} that
        allow \code{this.path()} to be used within the site-wide startup
        profile file or a user profile.}

      \item{\code{this.path(contents = TRUE)} in \sQuote{Rgui} will now return
        the contents of an unsaved document, or return the unsaved changes of a
        saved document.}

      \item{Added \code{set.gui.path()} for implementing \code{sys.path()} for
        arbitrary GUIs.}
    }
  }
}
\section{CHANGES IN VERSION 2.1.0 (2023-09-12)}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item{\code{this.path()} is now compatible with function
        \code{shiny::runApp()}.}

      \item{\code{this.path()} is now compatible with function
        \code{plumber::plumb()}.}

      \item{Added \code{tryCatch3()} which accepts condition expressions
        instead of condition handlers. Added active binding
        \codeorsQuote{last.condition} which will be the last thrown-and-caught
        condition in \code{tryCatch3()}.}

      \item{When \code{tryCatch2()} and \code{tryCatch3()} catch a condition
        but the respective handler / / expression is missing, the next
        non-missing handler / / expression is used instead.}
    }
  }
  \subsection{BUG FIXES}{
    \itemize{
      \item{\code{set.env.path()} and \code{set.src.path()} now return their
        input instead of returning \code{NULL}.}

      \item{\code{set.src.path()} now accepts an expression vector containing a
        source file.}
    }
  }
}
\section{CHANGES IN VERSION 2.0.0 (2023-08-08)}{
  \subsection{SIGNIFICANT USER-VISIBLE CHANGES}{
    \itemize{
      \item{\code{this.path()} was misleading users into believing it returns
        the path of the script in which it is written. In actuality, it returns
        the path of the executing script. As such, it was renamed to
        \code{sys.path()} in accordance with the other call stack inspection
        functions \code{sys.call()}, \code{sys.frame()}, \code{sys.nframe()},
        and \code{sys.function()}. Related functions that were renamed:
        \tabular{rll}{
               \tab old                   \tab new                  \cr
          [1,] \tab this.dir              \tab sys.dir              \cr
          [2,] \tab here                  \tab sys.here             \cr
          [3,] \tab this.proj             \tab sys.proj             \cr
          [4,] \tab rel2here              \tab rel2sys.dir          \cr
          [5,] \tab LINENO                \tab sys.LINENO           \cr
          [6,] \tab try.this.path         \tab try.sys.path         \cr
          [7,] \tab set.this.path         \tab set.sys.path         \cr
          [8,] \tab unset.this.path       \tab unset.sys.path       \cr
          [9,] \tab set.this.path.jupyter \tab set.sys.path.jupyter
        }

        \code{this.path()} was replaced by a function that returns the path of
        the script in which it is written. It does this by looking for a source
        reference with an appropriate source file, or inspects the environment
        in which it is called for an associated path, or returns the path of
        the executing script, as before. Related functions that were replaced:
        \tabular{rllll}{
          [1] \tab this.dir  \tab here       \tab ici           \tab this.proj  \cr
          [5] \tab rel2here  \tab LINENO     \tab try.this.path \tab check.path \cr
          [9] \tab check.dir \tab check.proj \tab               \tab
        }}
    }
  }
  \subsection{NEW FEATURES}{
    \itemize{
      \item{Added \code{env.path()} for determining the path associated with
        the top level environment. Also added \code{env.dir()},
        \code{env.here()}, \code{env.proj()}, \code{rel2env.dir()},
        \code{env.LINENO()}, and \code{try.env.path()}.}

      \item{Added \code{src.path()} for determining the path associated with
        its source reference. Also added \code{src.dir()}, \code{src.here()},
        \code{src.proj()}, \code{rel2src.dir()}, \code{src.LINENO()}, and
        \code{try.src.path()}.}

      \item{Added \code{set.env.path()} and \code{set.src.path()} to work along
        side \code{set.sys.path()}.}

      \item{\code{this.path()} now accepts arguments \codeorsQuote{n},
        \codeorsQuote{envir}, \codeorsQuote{matchThisEnv}, and
        \codeorsQuote{srcfile} to finely control how the path will be
        retrieved.}

      \item{Added \code{rel2proj()} for turning absolute paths into relative
        paths against the script's project root. Also added
        \code{rel2sys.proj()}, \code{rel2env.proj()}, and
        \code{rel2src.proj()}.}

      \item{\code{set.sys.path()} now accepts argument \codeorsQuote{ofile}
        specifying the original file argument. This overwrites the value
        returned by \code{sys.path(original = TRUE)}.}

      \item{Added active bindings \codeorsQuote{FILE} and \codeorsQuote{LINE}
        linking to \code{try.this.path()} and \code{LINENO()} to be used in a
        similar manner to the macros \codeorsQuote{__FILE__} and
        \codeorsQuote{__LINE__} in C.}

      \item{Added \code{path.functions()} which accepts a file path and
        constructs a set of path-related functions, similar to
        \code{this.path()} and associated.}
    }
  }
  \subsection{DEPRECATED AND DEFUNCT}{
    \itemize{
      \item{\code{inside.source()} and \code{set.this.path()} are defunct and
        should be replaced with \code{set.sys.path()}.}

      \item{\code{unset.this.path()} is defunct and should be replaced with
        \code{unset.sys.path()}.}

      \item{\code{set.this.path.jupyter()} is defunct and should be replaced
        with \code{set.sys.path.jupyter()}.}

      \item{\code{reset.this.proj()} is defunct and should be replaced with
        \code{reset.proj()}.}
    }
  }
  \subsection{BUG FIXES}{
    \itemize{
      \item{Fixed display issues in LaTeX and plain text documentation.}

      \item{\code{reset.proj()} now has environment
        \code{<environment: namespace:this.path>} as do all other exported
        functions.}

      \item{Fixed \code{sys.path(contents = TRUE)} in \sQuote{\RStudio} and
        \sQuote{\VSCode} having a trailing blank string.}

      \item{Fixed \code{sys.path(contents = TRUE)} in \sQuote{\Jupyter} if the
        executing script had yet to be determined.}

      \item{Fixed \code{sys.path(original = NA)} in \sQuote{\RStudio},
        \sQuote{\VSCode}, and \sQuote{\Jupyter} saying
        \code{missing value where TRUE/FALSE needed}.}

      \item{Fixed \code{relpath()} on Windows when comparing paths with network
        shares.}

      \item{Fixed \pkg{package:this.path} used in the site-wide startup profile
        file, a user profile, or the function \code{.First()} on Windows in
        \sQuote{\RStudio}.}
    }
  }
}
