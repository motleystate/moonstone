import logging
from unittest import TestCase

from moonstone.utils.log_messages import warn_once, reset_warnings_decorator

logger = logging.getLogger(__name__)


class TestLogMessage(TestCase):

    def test_warn_once_simple_case(self):
        with self.assertLogs('tests.utils.test_log_messages', level='WARNING') as log:
            warn_once(logger, "message to only print once.")
            warn_once(logger, "message to only print once.")
            self.assertEqual(len(log.output), 1)
            self.assertIn("WARNING:tests.utils.test_log_messages:message to only print once.", log.output)

    def test_warn_once_more_messages(self):
        def testing_function():
            warn_once(logger, "tracked warning message #1")
            warn_once(logger, "tracked warning message #1")
            logger.warning("other warning message")   # not tracked
            logger.warning("not tracked")             # not tracked
            for i in range(2, 11):
                warn_once(logger, f"tracked warning message #{i}")
            warn_once(logger, "tracked warning message #1")  # not printed because in the 10 messages tracked
            # But kick itself out to put itself back as most recent
            # -> I don't like it but I should create a personalized tracker and I don't care enough

            warn_once(logger, "other warning message")            # wasn't tracked so printed againt, but tracked now
            warn_once(logger, "tracked warning message #2")       # printed because not in the storage anymore

        with self.assertLogs('tests.utils.test_log_messages', level='WARNING') as log:
            testing_function()
            self.assertEqual(len(log.output), 14)
            self.assertEqual(log.output[0], "WARNING:tests.utils.test_log_messages:tracked warning message #1")
            self.assertEqual(log.output[1], "WARNING:tests.utils.test_log_messages:other warning message")
            self.assertEqual(log.output[2], "WARNING:tests.utils.test_log_messages:not tracked")
            self.assertEqual(log.output[3], "WARNING:tests.utils.test_log_messages:tracked warning message #2")
            # ...
            self.assertEqual(log.output[11], "WARNING:tests.utils.test_log_messages:tracked warning message #10")
            self.assertEqual(log.output[12], "WARNING:tests.utils.test_log_messages:other warning message")
            self.assertEqual(log.output[13], "WARNING:tests.utils.test_log_messages:tracked warning message #2")

    def test_reset_warnings_decorator(self):
        def testing_function_without_decorator(n):
            warn_once(logger, "tracked warning message (w/o decorator)")
            warn_once(logger, "tracked warning message (w/o decorator)")
            warn_once(logger, f"different log message {n}")   # because otherwise during the second calling I get
            # AssertionError: no logs of level WARNING or higher triggered on tests.utils.test_log_messages

        with self.assertLogs('tests.utils.test_log_messages', level='WARNING') as log:
            testing_function_without_decorator(1)
            self.assertEqual(len(log.output), 2)
            self.assertEqual(
                log.output[0], "WARNING:tests.utils.test_log_messages:tracked warning message (w/o decorator)")
            self.assertEqual(log.output[1], "WARNING:tests.utils.test_log_messages:different log message 1")

        with self.assertLogs('tests.utils.test_log_messages', level='WARNING') as log:
            testing_function_without_decorator(2)
            self.assertEqual(len(log.output), 1)     # cache not cleared so it doesn't print the first message
            # but it print the second message that is unique
            self.assertEqual(log.output[0], "WARNING:tests.utils.test_log_messages:different log message 2")

        @reset_warnings_decorator
        def testing_function_with_decorator():
            warn_once(logger, "tracked warning message (w/ decorator)")
            warn_once(logger, "tracked warning message (w/ decorator)")

        with self.assertLogs('tests.utils.test_log_messages', level='WARNING') as log:
            testing_function_with_decorator()
            self.assertEqual(len(log.output), 1)
            self.assertEqual(
                log.output[0], "WARNING:tests.utils.test_log_messages:tracked warning message (w/ decorator)")

        with self.assertLogs('tests.utils.test_log_messages', level='WARNING') as log:
            testing_function_with_decorator()
            self.assertEqual(len(log.output), 1)
            self.assertEqual(
                log.output[0], "WARNING:tests.utils.test_log_messages:tracked warning message (w/ decorator)")
