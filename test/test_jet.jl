using WaveguideQED, JET

using JET: ReportPass, BasicPass, InferenceErrorReport, UncaughtExceptionReport

# Custom report pass that ignores `UncaughtExceptionReport`
# Too coarse currently, but it serves to ignore the various
# "may throw" messages for runtime errors we raise on purpose
# (mostly on malformed user input)
struct MayThrowIsOk <: ReportPass end

# ignores `UncaughtExceptionReport` analyzed by `JETAnalyzer`
(::MayThrowIsOk)(::Type{UncaughtExceptionReport}, @nospecialize(_...)) = return

# forward to `BasicPass` for everything else
function (::MayThrowIsOk)(report_type::Type{<:InferenceErrorReport}, @nospecialize(args...))
    BasicPass()(report_type, args...)
end

rep = report_package("WaveguideQED";
    report_pass=MayThrowIsOk(), # TODO have something more fine grained than a generic "do not care about thrown errors"
    ignored_modules=(
        ))

@show rep
@test_broken length(JET.get_reports(rep)) == 0 # TODO JET reports issues
