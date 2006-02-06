def encode(data):
  edata = ""
  for i in xrange(len(data)):
    edata += "%.2x" % ord(data[i])
  return edata

def create_encoded(image_file_name):
  import wx
  img = wx.Image(name=image_file_name)
  w,h = img.GetSize()
  print 'img = img_data(width=%d, height=%d, mask=-1, encoded_data = """\\' % (
    w, h)
  encoded = encode(img.GetData())
  while (len(encoded) > 0):
    print encoded[:78]+"\\"
    encoded = encoded[78:]
  print '""")'

class img_data:

  def __init__(self, width, height, mask, encoded_data):
    self.width = width
    self.height = height
    self.data = self.decode(encoded_data)
    self.mask = mask * 3

  def get_width(self): return self.width
  def get_height(self): return self.height
  def get_size(self): return (self.width, self.height)
  def get_data(self): return self.data
  def get_mask(self): return self.mask

  def decode(self, edata):
    hex_chars = {"0":  0, "1":  1, "2":  2, "3":  3,
                 "4":  4, "5":  5, "6":  6, "7":  7,
                 "8":  8, "9":  9, "a": 10, "b": 11,
                 "c": 12, "d": 13, "e": 14, "f": 15}
    data = ""
    for i in xrange(0, len(edata), 2):
      data += chr(hex_chars[edata[i]] * 16 + hex_chars[edata[i+1]])
    return data

  def as_wx_Bitmap(self):
    import wx
    w,h = self.get_size()
    data = self.get_data()
    mask = self.get_mask()
    img = wx.EmptyImage(w, h)
    img.SetData(data)
    if (mask >= 0):
      img.SetMaskColour(ord(data[mask]), ord(data[mask+1]), ord(data[mask+2]))
      img.SetMask()
    return img.ConvertToBitmap()

center_img = img_data(width=20, height=20, mask=-1, encoded_data = """\
ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff\
ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff\
dadaff8b8bff5e5eff4d4dff5c5cff8a8affd9d9fffffffffffffffffffffffffff7f7f7ffffff\
ffffffffffffffffffffffffffffffe7e7ff5959ff0d0dff5454ff8c8cff9e9eff8c8cff5656ff\
0f0fff5454ffe4e4ffffffff9f9f9f2f2f2fffffffffffffffffffffffffffffffd0d0ff1a1aff\
6262ffe9e9fffffffffffffffffffffffffffdfdfdebebff6666ff1818fe6e6e87161616c6c6c6\
ffffffffffffffffffffffffe5e5ff1a1aff8f8fffffffffffffffffffffffffffffffffffffff\
9f9f9f999999fcfcfc3e3e6f000021bbbbd7ffffffffffffffffffffffffffffff5858ff6363ff\
ffffffffffffffffffffffffffffffffffffffffff3939390404043f3f3f2f2f2f5b5be45151ff\
ffffffffffffffffffffffffdadaff0e0effeaeaffffffffffffffffffffffffffffffffffffff\
d3d3d3010101000000000000b4b4b4ededff1010ffd4d4ffffffffffffffffffff8989ff5656ff\
ffffffffffffffffffffffffffffffffc7baff91756f524a000000000000000000131313dadada\
5c5cff8383ffffffffffffffffffff5a5aff8d8dffffffffffffffffffffffffffffc7b9ff3402\
f631001404000a09094545458b8b8bd1d1d1fdfdfd9595ff5555ffffffffffffffffffff4b4bff\
a1a1ffffffffffffffffffffffffffff8f73ff3300e22d00b72500f4896effffffffffffffffff\
ffffffa8a8ff4545ffffffffffffffffffff5a5aff9090ffffffffffffffffffffffffffffc3b4\
ff3401ff3300ff3401ffbfb0ffffffffffffffffffffffff9595ff5353ffffffffffffffffffff\
8787ff5959ffffffffffffffffffffffffffffffffffc1b2ff8a6cffbfafffffffffffffffffff\
ffffffffffff5f5fff8181ffffffffffffffffffffd5d5ff0e0effededffffffffffffffffffff\
fffffffffffffffffffffffffffffffffffffffffffffffff1f1ff1212ffd0d0ffffffffffffff\
ffffffffffff5151ff6b6bffffffffffffffffffffffffffffffffffffffffffffffffffffffff\
ffffffffffff7070ff4b4bffffffffffffffffffffffffffffffffe2e2ff1616ff9898ffffffff\
ffffffffffffffffffffffffffffffffffffffffffffffff9d9dff1515ffdedeffffffffffffff\
ffffffffffffffffffffffffcacaff1616ff6d6dffefefffffffffffffffffffffffffffffffff\
f1f1ff7272ff1414ffc6c6ffffffffffffffffffffffffffffffffffffffffffffffffffe0e0ff\
4e4eff1111ff5f5fff9696ffaaaaff9898ff6161ff1111ff4a4affdedeffffffffffffffffffff\
ffffffffffffffffffffffffffffffffffffffffffffffffd1d1ff8080ff5151ff4444ff5252ff\
7f7fffceceffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff\
ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff\
ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff\
ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff\
""")

spiral_img = img_data(width=20, height=20, mask=-1, encoded_data = """\
ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff\
fffffffffffffffffffffffffffefffffdfeffffffffffffffffffffffffffffffffffffffffff\
ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff\
fffcfdffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff\
fffffffffffffffffffffffffffffffffffffffffefcf8f9fffffffffffffefffcfefffeffffff\
fffffffdfffffefefffcfffffdfffffffdfffffdfffffefffefffffefffffefffffdfffffeffff\
fefffefdfffefffffffefdfbfefcfcfefefdfbfffdfffefdfffefdfffefefdfffffefffffffeff\
fffefffffefffdfffffdfffffffefffdfffffefffffbfffffefefffffffffffefefffdfffdffff\
fefefefefefdfcfefefefafefaeef3eac2ceda9aadd7839ad58a9fe0a9baf5dbe5fcfafbfefdfc\
fcfefefffffffefdfffffffffffffffefffefdfffffefcfffefafbf3d7decc738db5375cc14d6e\
c77187ce7891c96381b62a51b41640c5647af7ebedfefefefffffffdfefffffffffffffffdfefb\
fefdfdfdf6fad593a7ad2c50d4899ff8e4e9fdf6f7f7e5eaf4dfe5feeef5fbe9f2ca6886b10a39\
c96e88fff9fbfffffffffffffffffffffffffdfdfdfef5f7cf6f8bb11a47eabac9fbfcfce7b5c7\
bb3e66b31747b2173dba274ec2607df5ccd8af143ebc4968fdfbfafffffffefffeffffffffffff\
fffbfed994a7b40b38d38097fcfbfcf6dee7ab1f4abd345bc44b6dbe3c5fc95e7ff2d6e0e1a1b3\
b20c39d8899fe5c3cefffffffefcfdfffffffffffffaebf1b1284fb1103cefc8d6fdfefef0d2d9\
ab2249f9d8e1fef9fbfef6f8fdf0f3dfa4b3ae1f46c7607ae6afc0cc89a1fffffffffcffffffff\
ffffffedc2cdb00c39af123ceecfd9fcfdfdfcf8f7cd738cb34361cb768bcc768cb84363b83e61\
dea0b2e7b4c3b6546ffae9eeffffffffffffffffffffffffebb9c7af0c38b10a36d07b93fcf8fa\
fefefefcfafbf3dbe3e2b9c4e3b3c2f1d4dbf8e0e4d17d94b23a5bf3d8dffffffeffffffffffff\
fffffffffffff8e1eab01f48b50a35b0113acf7b90f3d2dafeecf3feedf4f6dfe5e8bdc4d0778e\
ae2a50ba4e6ff5dce4fffffffffffffffffffffffffffffffffffffefcfde0a3b6af1f45b20a35\
b40a36b1123bb32048ae2249ae1841ae143db9385edfa1b6fdf5f9fffeffffffffffffffffffff\
fffffffffffffffffffbfefefafcfcf2d7ded5889bc35874be4666bf4b69c8657ed98c9eecc3d0\
fcf5f7fdfefdfffffffffffffffffffffffffffffffffffffffffffffffffffffffbfefffefdff\
fefefefffefefffdfefffbfdfffcfdfffdfefefefefefbfdfffdfdffffffffffffffffffffffff\
fffffffffffffffffffffffffffffffefefffcffffffffffffffffffffffffffffffffffffffff\
fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffefd\
fefffdffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff\
ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff\
ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff\
""")

fit_img = img_data(width=20, height=20, mask=0, encoded_data = """\
ece9d8ece9d8ece9d8ece9d8ece9d8ece9d8ece9d8ece9d8ece9d8ece9d8ece9d8ece9d8ece9d8\
ece9d8ece9d8ece9d8ece9d8ece9d8ece9d8ece9d8ece9d8ece9d8ece9d8ece9d8ece9d8ece9d8\
ece9d8ece9d8ece9d8ece9d8ece9d8ece9d8ece9d8ece9d8ece9d8ece9d8ece9d8ece9d8ece9d8\
ece9d8528bc56694c992a9d096add394acd498afd59ab2d79db3d7a0b7d9a5bddfa8c2e5aac3e5\
adc5e6aec7e7b2cae8b7cdeabfd2ecbed3eca8c7e798bde46191c4e3edf7e3edf7e3edf7e3edf7\
e3edf7e3edf7e3edf7c1d4e734679ac1d4e7e3edf7e3edf7e3edf7e3edf7e3edf7e3edf7e3edf7\
e3edf798bde46791c1e1ebf6e4edf8e7eff9ebf2faeef3faf1f6fcc1d4e734679a34679a34679a\
c1d4e7ffffffffffffffffffffffffffffffffffffffffff98bde3688dbcdeeaf6e2ebf7e5eef8\
e9f0f9ecf3faeff5fb34679a34679a34679a34679a34679afefeffffffffffffffffffffffffff\
ffffffffffff98bde36a8ab8dbe8f6e0eaf7e3ecf8e7eff8eaf1f9eef3faf0f5fbf3f7fc34679a\
f8fbfdfafcfefcfdfeffffffffffffffffffffffffffffffffffff97bde36b8ab7d9e6f5dde9f6\
e0ebf7e4edf8e8eff8ebf1faeef4faf1f6fc34679af7fafdf9fbfdfcfcfefdfefeffffffffffff\
ffffffffffffffffff9abee36b8ab7d7e5f5c1d4e734679adde9f5e6edf8e9f0f9ecf2faeff5fa\
34679af5f8fcf8fbfdfafcfdfcfdfefefeffb9cbdc4674a3ffffffffffffa8c3e46b8ab7b3c9e1\
36689b3c6d9edfe9f7e3ecf8e7eff8eaf1faedf4fbf1f5fbf3f7fcf6f9fcf8fbfdfbfcfefdfefe\
c6d4e33366994d7aa6ffffffaac3e46b8ab734679a33669933669934679a34679a35679ae8f0f9\
ebf2f9eff4fbf2f6fbf4f8fdf7f9fd35679a35679a35679a33669933669935679aaac4e56988b6\
a0bbd73366993b6c9ec8d9ebd0deeedfeaf6e5eef8e9f1f9a0bbd7f0f4fbf2f6fcf5f8fceff4f9\
eaf0f6b6c8db3366993d6d9eeef2f7aac5e56d8bb9cbdcf1aac2dd336699d6e4f4dce7f6e0eaf6\
e3ecf8e6eef834679aedf3faf0f5fbf3f7fcf6f9fdf9fbfeb4c7da37699bebf1f5ffffffaac5e5\
7795bfccddf1cedef2d2e1f3d6e3f4d9e6f5dde8f6e1ebf7e4edf734679aebf2f9eff4faf2f6fb\
f4f7fcf6fafdf9fbfdfbfcfefdfefeffffffabc5e57895c0c9dbf0ccddf1d0e0f3d3e2f4d7e4f4\
dae7f534679a34679a34679a34679a34679aeff4fbf2f7fcf5f9fcf8fafdfafcfdfdfdfeffffff\
abc5e5718ebbc7daf0cadcf1cddef2d1e0f3d4e3f4d8e5f4c1d4e734679a34679a34679ac1d4e7\
edf3fbf0f5fbf3f7fcf6f9fdf9fbfdfafcfeffffffabc6e56b8ab9c5d9f0c8dbf0cbdcf1cfdff2\
d2e1f3d6e3f4dae6f5c1d4e734679ac1d4e7e7eff9ebf2faeef3faf1f6fbf4f8fcf6f9fdf9fbfd\
ffffffacc7e76887b75f82b66f8ebd7393bf7595c17896c37b99c67f9dc984a1cc87a5d08baad2\
8facd694b1d897b4db9cb8dda0bce0a4c0e2a8c3e4aac5e6aac5e6ece9d8ece9d8ece9d8ece9d8\
ece9d8ece9d8ece9d8ece9d8ece9d8ece9d8ece9d8ece9d8ece9d8ece9d8ece9d8ece9d8ece9d8\
ece9d8ece9d8ece9d8ece9d8ece9d8ece9d8ece9d8ece9d8ece9d8ece9d8ece9d8ece9d8ece9d8\
ece9d8ece9d8ece9d8ece9d8ece9d8ece9d8ece9d8ece9d8ece9d8ece9d8\
""")

align_img = img_data(width=20, height=20, mask=-1, encoded_data = """\
cfffffccffffc6ffffc1ffffbdffffb7ffffb2ffffadffffa8ffffa3ffff9dffff98ffff93ffff\
8effff88ffff84ffff7effff79ffff74ffff6fffffccffff0886ff0886ffbcffffb7ffffb2ffff\
adffffa7ffffa2ffff9dffff98ffff93ffff8dffff89ffff84ffff7effff79ffff74ffff6fffff\
69ffffc7ffff0886ff0886ffb7ffffff0000ff0000ff0000ff0000ff0000ff0000ff0000ff0000\
ff0000ff0000ff0000ff0000ff0000ff000069ffff63ffffc2ffff0886ff0886ffb2ffffff0000\
a8ffffa2ffff9dffff98ffff93ffff8dffff88ffff83ffff7effff79ffff74ffff6fffffff0000\
64ffff5effffbcffff0886ff0886ffacffffff0000a1ffff9dffff97ffff93ffff8effff88ffff\
82ffff7effff79ffff73ffff6fffff69ffffff00005effff59ffffb7ffff0886ff0886ffa7ffff\
ff00009dffff97ffff93ffff8dffff88ffff83ffff7effff78ffff73ffff6effff69ffff64ffff\
ff000059ffff53ffffb1ffff0886ff0886ffa1ffffff000097ffff92ffff8dffff88ffff83ffff\
7effff78ffff73ffff6effff68ffff63ffff5effffff000054ffff4effffacffff0886ff0886ff\
9cffffff000092ffff8dffff88ffff82ffff7dffff78ffff73ffff6dffff69ffff63ffff5effff\
58ffffff00004effff48ffffa7ffff0886ff0886ff97ffffff00008dffff87ffff83ffff7dffff\
78ffff73ffff6effff68ffff63ffff5dffff58ffff53ffffff000048ffff43ffffa1ffff0886ff\
0886ff92ffffff000087ffff82ffff7dffff78ffff73ffff6dffff68ffff63ffff5dffff58ffff\
52ffff4dffffff000042ffff3dffff9cffff0886ff0886ff8cffffff000082ffff7dffff78ffff\
72ffff6effff68ffff63ffff5effff58ffff53ffff4dffff48ffffff00003cffff37ffff97ffff\
0886ff0886ff87ffffff00007dffff77ffff72ffff6dffff68ffff62ffff5dffff58ffff52ffff\
4dffff48ffff42ffffff000038ffff32ffff91ffff0886ff0886ff82ffffff000078ffff72ffff\
6dffff68ffff62ffff5dffff58ffff52ffff4dffff47ffff42ffff3dffffff000032ffff2dffff\
8cffff0886ff0886ff7cffffff000072ffff6cffff67ffff63ffff5cffff58ffff52ffff4cffff\
47ffff42ffff3dffff37ffffff00002cffff27ffff86ffff0886ff0886ff77ffffff00006cffff\
68ffff62ffff5dffff58ffff52ffff4cffff47ffff41ffff3cffff37ffff31ffffff000027ffff\
21ffff81ffff0886ff0886ff72ffffff0000ff0000ff0000ff0000ff0000ff0000ff0000ff0000\
ff0000ff0000ff0000ff0000ff0000ff000021ffff1bffff7cffff0886ff0886ff6cffff68ffff\
62ffff5dffff57ffff51ffff4dffff46ffff42ffff3cffff36ffff31ffff2bffff26ffff21ffff\
1cffff16ffff77ffff0886ff0886ff0886ff0886ff0886ff0886ff0886ff0886ff0886ff0886ff\
0886ff0886ff0886ff0886ff0886ff0886ff0886ff0886ff11ffff72ffff0886ff0886ff0886ff\
0886ff0886ff0886ff0886ff0886ff0886ff0886ff0886ff0886ff0886ff0886ff0886ff0886ff\
0886ff0886ff0bffff6cffff67ffff62ffff5cffff57ffff51ffff4bffff46ffff41ffff3bffff\
36ffff30ffff2cffff26ffff20ffff1cffff16ffff10ffff0bffff05ffff\
""")
