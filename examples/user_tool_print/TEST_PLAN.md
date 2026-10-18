# Print Import acceptance test

User input: original.py, a single print statement.
AI adaptation: adapter.py wraps that statement in the v1 run entrypoint.
Expected complete feature: import -> user tool registered -> execution action ->
stdout ASTRAFEX_USER_PRINT_OK visible in run log, on browser and standalone.
A static_pass_needs_review report alone does not meet these acceptance criteria.
