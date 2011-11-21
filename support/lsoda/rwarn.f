      subroutine rwarn(msg)
      character*(*) msg
        call rwarnc(msg, len(msg))
      end

