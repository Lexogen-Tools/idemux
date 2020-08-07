import idemux.idemux as demux
import cProfile
cProfile.run('demux.main()', sort='cumulative')


