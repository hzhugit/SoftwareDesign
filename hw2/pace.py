# If I leave my house at 6:52 am and run 1 mile at an easy pace (8:15 per mile), then 3 miles at tempo (7:12 per mile) and 1 mile at easy pace again, what time do I get home for breakfast?

import math
import unum

start_time = 6*60 + 52
easy_pace = 8.0 + 15/60.0
tempo = 7 + 12/60.0

total_time = 2 * easy_pace + 3 * tempo
final_time = start_time + total_time
print divmod(final_time, 60)
