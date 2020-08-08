import idemux.idemux as demux
import cProfile
#demux.main()
cProfile.run('demux.main()', sort='cumulative')


