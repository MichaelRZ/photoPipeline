#include <iostream>
#include <vector>
#include <string>

struct{
  vec3 radiance;
  float NV;
  uint photoColor;
}texelStruct;

int main(){
  std::vector<float> minDepths; //index x + y*frame_width, {screenspace{depth}}
  std::vector<std::vector<texelStruct>> texelDatas; // {uv space {camera angles {texel solving info}}}
  std::vector<std::vector<uint8_t>> result; // {uv space {texture channels}}
  for each camera angle{
    depthPass(); // fills minDepths
    renderPasses(); // fills texelDatas, empty minDepths
  }
  solves(); //fills result, empty texelDatas
  //write to texture files here

  return 0;

}