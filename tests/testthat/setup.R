#Open a null graphics device for the whole test run.
#
#Several code paths need a device even though they never draw for the user:
#`love.plot()` with more than one statistic converts each plot to a grob to
#assemble the gtable, and grob conversion needs a device to measure text. With no
#device open R starts the default one, which writes an `Rplots.pdf` into
#tests/testthat/. Opening the null device here means nothing is ever written.
#
#No teardown is registered: the device closes when the test process exits, and
#tests that need their own device use `local_null_device()` (see helpers.R), which
#opens and closes one around itself.
grDevices::pdf(NULL)
