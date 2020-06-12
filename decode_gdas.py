import sys
import ncepbufr
import numpy as np

#max_subs_per_msg = 1
max_subs_per_msg = sys.maxsize 

#max_real_msgs = 10
max_real_msgs = sys.maxsize

bufr = ncepbufr.open(sys.argv[1])
real_msg_count = 0
while bufr.advance() == 0:
    print('msg_type: ', bufr.msg_type,
          '- msg_date: ', bufr.msg_date,
          '- msg_counter: ', bufr.msg_counter,
          '- receipt_time: ', bufr.receipt_time,
          '- subsets: ', bufr.subsets, file=sys.stderr)
    if bufr.subsets > 0 and real_msg_count < max_real_msgs:
        subs_count = 0
        while bufr.load_subset() == 0:
            print('real_msg_count =', real_msg_count, ', subs_count =',
                  subs_count, file=sys.stderr)
            bufr.print_subset()
            subs_count += 1
            if subs_count > max_subs_per_msg:
                break
    if bufr.subsets > 0:
        real_msg_count += 1
    if real_msg_count >= max_real_msgs:
        break

bufr.close()
