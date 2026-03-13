# Adding tool to galaxy

## Add tool definition
Modify the tool_conf.xml file to include the new tool. This file is typically located in the Galaxy installation directory under config/tool_conf.xml or a similar path.

```xml
<!-- lib/galaxy/config/sample/tool_conf.xml.sample -->
  <section id="uont" name="uONT">
    <tool file="uont/uont_assemble.xml" />
  </section>
```

## Add the definition to the tools directory

Make the `tools/uont` directory if it doesn't exist and copy the `uont_assemble.xml` file to the appropriate location in the Galaxy tools directory. This file defines the interface for the uONT assembly tool in Galaxy.